#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>

#include "bamHelpers.hpp"
#include "version.hpp"

#ifndef kroundup32
#define kroundup32(x) (--(x), (x) |= (x) >> 1, (x) |= (x) >> 2, (x) |= (x) >> 4, (x) |= (x) >> 8, (x) |= (x) >> 16, ++(x))
#endif

// vector2cigar will convert a vectorised CIGAR back to a htslib compatible one.
uint32_t* vector2cigar(std::vector<uint32_t>* cigarVector, uint32_t cigarLen)
{
    uint32_t* newCigar_ptr = (uint32_t*)malloc(sizeof(uint32_t) * cigarLen);
    uint32_t* newCigar = newCigar_ptr;
    for (auto it = cigarVector->begin(); it != cigarVector->end(); ++it)
    {
        // collect both the CIGAR operation and its length from successive vector slots
        uint32_t cigOp = *it;
        ++it;
        uint32_t cigLen = *it;

        // pack both the operation and length into a single uint32
        cigLen = cigLen << BAM_CIGAR_SHIFT;
        *newCigar_ptr = (cigOp | cigLen);
        newCigar_ptr++;
    }
    return newCigar;
}

// replaceCigar will replace a records current CIGAR with a new one.
void replaceCigar(bam1_t* record, uint32_t* cigar, uint32_t length)
{
    // if the length of the new matches the old, we don't need to mess with memory re-allocations
    if (length == record->core.n_cigar)
    {
        memcpy(record->data + record->core.l_qname, cigar, length * 4);
        return;
    }

    // otherwise, reallocate the record data if the new CIGAR will overflow the current maximum length of bam1_t::data
    if (record->l_data + (length - record->core.n_cigar) * 4 > record->m_data)
    {
        record->m_data = record->l_data + (length - record->core.n_cigar) * 4;
        kroundup32(record->m_data);
        record->data = (uint8_t*)realloc(record->data, record->m_data);
    }

    // out with the old and in with the new
    int oldGuff = record->core.l_qname + record->core.n_cigar * 4;
    memmove(record->data + record->core.l_qname + length * 4, record->data + oldGuff, record->l_data - oldGuff);
    memcpy(record->data + record->core.l_qname, cigar, length * 4);

    // update the length of data and length of CIGAR for the record
    record->l_data += (length - record->core.n_cigar) * 4;
    record->core.n_cigar = length;
}

// reverseCigar reverses a CIGAR.
// use XOR to swap bits without using an intermediate variable
// see: https://www.geeksforgeeks.org/swap-two-numbers-without-using-temporary-variable/
void reverseCigar(uint32_t* cigar, int length)
{
    for (int i = 0; i < length / 2; ++i)
    {
        cigar[i] ^= cigar[(length - 1) - i];
        cigar[(length - 1) - i] ^= cigar[i];
        cigar[i] ^= cigar[(length - 1) - i];
    }
}

// primerTrim will perform the softmasking and return a new CIGAR (the caller must replace the existing CIGAR).
// this function essentially mirrors the python `trim(segment, primer_pos, end, debug)` func from align_trim.
void primerTrim(std::vector<uint32_t>* cigar, bam1_t* record, uint32_t maskEnd, bool reverseMask)
{

    // get the segment position in the reference (depends on if start or end of the segment is being processed)
    uint32_t pos;
    if (!reverseMask)
    {
        pos = record->core.pos;

        // if forward primer, use the start of the record and flip the CIGAR for easy appending
        std::reverse(cigar->begin(), cigar->end());
    }
    else
    {
        pos = bam_endpos(record);
    }

    // process the CIGAR to determine how much softmasking is required
    uint32_t eaten = 0;
    while (1)
    {
        uint32_t cigarOp;
        uint32_t cigarLen;

        // chomp CIGAR operations from the start/end of the CIGAR
        if (cigar->size() >= 2)
        {
            if (reverseMask)
            {
                cigarLen = cigar->back();
                cigar->pop_back();
                cigarOp = cigar->back();
                cigar->pop_back();
            }
            else
            {
                cigarOp = cigar->back();
                cigar->pop_back();
                cigarLen = cigar->back();
                cigar->pop_back();
            }
        }
        else
        {
            // ran out of cigar during soft masking - completely masked read will be ignored
            break;
        }

        // if the CIGAR operation consumes the reference sequence, increment/decrement the position by the CIGAR operation length
        if (bam_cigar_type(cigarOp) & 2)
        {
            (!reverseMask) ? pos += cigarLen : pos -= cigarLen;
        }

        // if the CIGAR operation consumes the query sequence, increment the number of CIGAR operations eaten by the CIGAR operation length
        if (bam_cigar_type(cigarOp) & 1)
            eaten += cigarLen;

        // stop processing the CIGAR if we've gone far enough to mask the primer
        if (!reverseMask && (pos >= maskEnd) && (cigarOp == BAM_CMATCH))
            break;
        if (reverseMask && (pos <= maskEnd) && (cigarOp == BAM_CMATCH))
            break;
    }

    // calculate how many extra matches are needed in the CIGAR
    uint32_t extra = (pos > maskEnd) ? pos - maskEnd : maskEnd - pos;
    if (extra)
    {
        if (reverseMask)
        {
            cigar->push_back(0);
            cigar->push_back(extra);
        }
        else
        {
            cigar->push_back(extra);
            cigar->push_back(0);
        }
        eaten -= extra;
    }

    // check we have smomething to softclip
    if (eaten <= 0)
        throw std::runtime_error("invalid cigar operation created - possibly due to INDEL in primer");

    // softmask the reverse/forward primer
    if (reverseMask)
    {
        cigar->push_back(BAM_CSOFT_CLIP);
        cigar->push_back(eaten);
    }
    else
    {
        // update the position of the leftmost mappinng base
        record->core.pos = pos - extra;

        // if proposed softmask leads straight into a deletion, shuffle leftmost mapping base along and ignore the deletion
        if (cigar->back() == BAM_CDEL)
        {
            while (1)
            {

                if (cigar->back() != BAM_CDEL)
                    break;

                // remove the deletion operation, add its length to the record start pos and then remove the length from the CIGAR as well
                cigar->pop_back();
                record->core.pos += cigar->back();
                cigar->pop_back();
            }
        }

        // add the soft clip
        cigar->push_back(eaten);
        cigar->push_back(BAM_CSOFT_CLIP);
    }

    // if we flipped the CIGAR at the start, flip it again
    if (!reverseMask)
        std::reverse(cigar->begin(), cigar->end());

    // check the the start/end operation of CIGAR has valid length
    if (cigar->at(2) <= 0 || cigar->end()[-2] <= 0)
        throw std::runtime_error("invalid cigar operation created - possibly due to INDEL in primer");
}

// TrimAlignment will softmask an alignment from its start/end up to the provided position.
void artic::TrimAlignment(bam1_t* record, unsigned int maskEnd, bool reverse)
{

    // copy the CIGAR to a vector
    uint32_t* origCigar = bam_get_cigar(record);
    std::vector<uint32_t> cigarVector;
    for (uint32_t i = 0; i < record->core.n_cigar; i++)
    {
        uint32_t cigOp = bam_cigar_op(origCigar[i]);
        uint32_t cigLen = bam_cigar_oplen(origCigar[i]);
        cigarVector.push_back(cigOp);
        cigarVector.push_back(cigLen);
    }

    // softmask the record
    primerTrim(&cigarVector, record, maskEnd, reverse);

    // convert updated CIGAR vector back to an array
    uint32_t cigarLen = cigarVector.size() / 2;
    auto newCigar = vector2cigar(&cigarVector, cigarLen);

    // replace the records original CIGAR with the updated one
    replaceCigar(record, newCigar, cigarLen);
}

// addLineToHeader is the function that actually updates the header and is called by the AddXXtoHeader functions.
void addLineToHeader(bam_hdr_t** header, const char* line)
{

    // get the required size for the new header
    size_t headerSize = std::strlen((*header)->text) + std::strlen(line) + 1;

    // copy the old header text into the new holder
    char* newHeaderText = (char*)malloc(headerSize);
    memcpy(newHeaderText, (*header)->text, std::strlen((*header)->text));

    // add a terminator to the old text and then add the new line
    newHeaderText[std::strlen((*header)->text)] = '\0';
    std::strcat(newHeaderText, line);

    // free the older header text and replace it with the new one
    free((*header)->text);
    (*header)->text = newHeaderText;
    newHeaderText = NULL;
    (*header)->l_text = headerSize - 1;
}

// findPreviousProg will check a BAM header for the last tool in the tool chain.
std::string findPreviousProg(bam_hdr_t** header)
{
    std::string PPtag;
    auto headerText = std::string((*header)->text);
    auto found = headerText.rfind("\tPN:");
    if (found != std::string::npos)
    {
        PPtag = "\tPP:";
        for (char const& c : headerText.substr(found + 3))
        {
            if (c == '\t' || c == ' ')
                break;
            PPtag.push_back(c);
        }
    }
    return PPtag;
}

// AddPGtoHeader will add a PG tag to an existing BAM header based on the provided command.
void artic::AddPGtoHeader(bam_hdr_t** header, const std::string& userCmd)
{
    std::string line = "@PG\tID:" + PROG_NAME + "\tPN:" + PROG_NAME;

    // see if we need to add a previous tool
    std::string pp = findPreviousProg(header);
    if (!pp.empty())
    {
        line.append(pp);
    }

    // add the prog version
    line.append("\tVN:" + GetVersion());

    // add the command
    line.append("\tCL:" + PROG_NAME + " " + userCmd + "\n\0");

    // update the header
    addLineToHeader(header, line.c_str());
}

// AddRGtoHeader will add a RG tag to an existing BAM header based on the provided read group.
void artic::AddRGtoHeader(bam_hdr_t** header, std::string& rg)
{
    std::string line = "@RG\tID:";
    line.append(rg + "\tPG:" + PROG_NAME + "\n\0");
    addLineToHeader(header, line.c_str());
}

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>

#include "bamHelpers.hpp"
#include "version.hpp"

#ifndef kroundup32
#define kroundup32(x) (--(x), (x) |= (x) >> 1, (x) |= (x) >> 2, (x) |= (x) >> 4, (x) |= (x) >> 8, (x) |= (x) >> 16, ++(x))
#endif

// replaceCIGAR will replace a record's current CIGAR with a new one.
void replaceCIGAR(bam1_t* record, uint32_t* cigar, uint32_t length)
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

// reverseCIGAR reverses a CIGAR.
// use XOR to swap bits without using an intermediate variable
// see: https://www.geeksforgeeks.org/swap-two-numbers-without-using-temporary-variable/
void reverseCIGAR(uint32_t* cigar, int length)
{
    for (int i = 0; i < length / 2; ++i)
    {
        cigar[i] ^= cigar[(length - 1) - i];
        cigar[(length - 1) - i] ^= cigar[i];
        cigar[i] ^= cigar[(length - 1) - i];
    }
}

// findQueryStart uses the CIGAR to find the mask end relative to the query sequence (instead of the ref sequence).
int32_t findQueryStart(uint32_t* cigar, uint32_t cigarLength, int32_t maskEnd, int32_t refStart)
{
    // set holders for the CIGAR operation and its length
    int cigOp;
    int32_t opLen;

    // set chomp counters for query and reference sequences
    int32_t qChomper = 0;
    int32_t rChomper = refStart;

    // iterate through the CIGAR and use bam_cigar_type to determine if reference or query consuming
    // see: https://github.com/samtools/htslib/blob/29d06c2f84dc20fbc5f94ca1b181292433b69a21/htslib/sam.h#L122
    for (uint32_t i = 0; i < cigarLength; ++i)
    {
        cigOp = bam_cigar_op(cigar[i]);
        opLen = bam_cigar_oplen(cigar[i]);

        // process reference consuming operation
        if (bam_cigar_type(cigOp) & 2)
        {

            // if teh reference has been chomped enough, we're ready to return the current query position
            if (maskEnd <= rChomper + opLen)
            {

                // if query consumed too, increment the query chomper before returning it
                if (bam_cigar_type(cigOp) & 1)
                    qChomper += (maskEnd - rChomper);
                return qChomper;
            }
            rChomper += opLen;
        }

        // otherwise we have a query consumer so chomp some more of the query
        if (bam_cigar_type(cigOp) & 1)
            qChomper += opLen;
    }
    return qChomper;
}

// primerTrim will perform the softmasking and return a new CIGAR (the caller must replace the existing CIGAR).
// TODO: will make this mirror the python `trim(segment, primer_pos, end, debug)` func from align_trim, just playing with variations for now.
artic::cigarHolder primerTrim(bam1_t* record, int32_t maskEnd, bool reverseMask)
{

    // get the current CIGAR and set up a holder for the new one (add an extra slot as CIGAR can could grow by 1 op type)
    uint32_t* cigar = bam_get_cigar(record);
    uint32_t* newCigar = (uint32_t*)malloc(sizeof(uint32_t) * (record->core.n_cigar + 1));

    // calc the max mask length possible
    int maxMaskLength = 0;
    (reverseMask) ? maxMaskLength = bam_cigar2qlen(record->core.n_cigar, bam_get_cigar(record)) - findQueryStart(cigar, record->core.n_cigar, maskEnd, record->core.pos) - 1
                  : maxMaskLength = findQueryStart(cigar, record->core.n_cigar, maskEnd, record->core.pos);

    // shouldn't happen, but allow for record only spanning primer region
    maxMaskLength = (maxMaskLength > 0) ? maxMaskLength : 0;

    // if we're masking the reverse, flip the CIGAR to make things easier
    reverseCIGAR(cigar, record->core.n_cigar);

    // set start position and reference trackers
    bool startPosReached = false;
    int32_t startPos = 0;
    int32_t refIncrement = 0;

    // set holders for the CIGAR operation and its length
    int cigOp;
    int32_t opLen;

    // set up two iterators for swapping between CIGARS
    uint32_t i = 0, j = 0;

    // process the CIGAR operations
    while (i < record->core.n_cigar)
    {

        // query bases have been masked, transfer over the CIGAR op and continue
        if (maxMaskLength == 0 && startPosReached)
        {
            newCigar[j] = cigar[i];
            i++;
            j++;
            continue;
        }

        // grab the op
        cigOp = bam_cigar_op(cigar[i]);
        opLen = bam_cigar_oplen(cigar[i]);

        // once masking is done, wait until we increment the startPos so it consumes both query and ref
        if (maxMaskLength == 0 && (bam_cigar_type(cigOp) & 1) && (bam_cigar_type(cigOp) & 2))
        {
            startPosReached = true;
            continue;
        }
        refIncrement = opLen;

        // consume query sequence
        if ((bam_cigar_type(cigOp) & 1))
        {

            // add the soft mask
            if (maxMaskLength >= opLen)
            {
                newCigar[j] = bam_cigar_gen(opLen, BAM_CSOFT_CLIP);
            }
            else if (maxMaskLength < opLen && maxMaskLength > 0)
            {
                newCigar[j] = bam_cigar_gen(maxMaskLength, BAM_CSOFT_CLIP);
            }
            else if (maxMaskLength == 0)
            {
                newCigar[j] = bam_cigar_gen(opLen, BAM_CSOFT_CLIP);
                j++;
                i++;
                continue;
            }
            j++;

            // update the mask progress
            refIncrement = std::min(maxMaskLength, opLen);
            auto oldLen = opLen;
            opLen = std::max(opLen - maxMaskLength, 0);
            maxMaskLength = std::max(maxMaskLength - oldLen, 0);
            if (opLen > 0)
            {
                newCigar[j] = bam_cigar_gen(opLen, cigOp);
                j++;
            }

            // check to see if we've masked enough
            if (maxMaskLength == 0 && (bam_cigar_type(newCigar[j - 1]) & 1) && (bam_cigar_type(newCigar[j - 1]) & 2))
                startPosReached = true;
        }

        // consume reference sequence
        if ((bam_cigar_type(cigOp) & 2))
            startPos += refIncrement;

        i++;
    }

    // if we were masking the reverse, we will have reversed the CIGAR to help, so flip it back again now
    if (reverseMask)
        reverseCIGAR(newCigar, j);

    // return a cigarHolder with the new CIGAR deets
    return {
        startPos,
        j,
        newCigar};
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

// TrimAlignment will softmask an alignment from its start/end up to the provided position.
void artic::TrimAlignment(bam1_t* record, unsigned int maskEnd, bool reverse)
{
    cigarHolder newCigar = primerTrim(record, maskEnd, reverse);
    record->core.pos += newCigar.startPos;
    replaceCIGAR(record, newCigar.cigar, newCigar.length);
    free(newCigar.cigar);
    return;
}

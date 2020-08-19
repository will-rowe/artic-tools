#include <fstream>
#include <htslib/sam.h>
#include <iostream>
#include <sstream>
#include <string>

#include "bamHelpers.hpp"
#include "log.hpp"
#include "softmask.hpp"

// getErrorMsg returns the error message for the softmasker error codes.
const char* getErrorMsg(int errorCode)
{
    switch (errorCode)
    {
        case 0:
            return "no error";
        case 1:
            return "softmasker is uninitialised";
        case 2:
            return "skipped as unmapped";
        case 3:
            return "skipped as supplementary";
        case 4:
            return "skipped as poor quality";
        default:
            return "unknown error";
    }
}

// _checkRecord returns an error if the currently held record fails filters and should be skipped.
unsigned int artic::Softmasker::_checkRecord(void)
{
    if (!_curRec)
        return 1;
    if (_curRec->core.flag & BAM_FUNMAP)
        return 2;
    if (_curRec->core.flag & BAM_FSUPPLEMENTARY)
        return 3;
    if (_curRec->core.qual < _minMAPQ)
        return 4;
    return 0;
}

// _getAmpliconCount returns the number of times the queried amplicon has been seen before.
unsigned int artic::Softmasker::_getAmpliconCount(void)
{
    std::string ampliconCounterKey = _curAmplicon->GetID();
    if (_curRec->core.flag & BAM_FREVERSE)
    {
        ampliconCounterKey.append("_reverse");
    }
    auto mapIter = _ampliconCounter.find(ampliconCounterKey);
    if (mapIter == _ampliconCounter.end())
    {
        _ampliconCounter.emplace(ampliconCounterKey, 1);
        return 1;
    }
    return ++mapIter->second;
}

// _reportLine will add the primer information for the current record to the open report.
// if verbose, the line will be logged to STDERR as well.
void artic::Softmasker::_reportLine(bool verbose)
{

    // calc Primer1Start and Primer2Start (for report compatability with Python code)
    auto maxSpan = _curAmplicon->GetMaxSpan();
    auto p1Start = std::abs(maxSpan.first - _curRec->core.pos);
    auto p2Start = std::abs(maxSpan.second - bam_endpos(_curRec));

    // get the info ready
    std::ostringstream buffer;
    buffer << bam_get_qname(_curRec) << "\t" << _curRec->core.pos << "\t" << bam_endpos(_curRec) << "\t" << _curAmplicon->GetID() << "\t" << _curAmplicon->GetForwardPrimer()->GetID() << "\t" << p1Start << "\t" << _curAmplicon->GetReversePrimer()->GetID() << "\t" << p2Start << "\t";
    (_curRec->core.flag & BAM_FSECONDARY) ? buffer << "True\t" : buffer << "False\t";
    (_curRec->core.flag & BAM_FSUPPLEMENTARY) ? buffer << "True\t" : buffer << "False\t";
    buffer << maxSpan.first << "\t" << maxSpan.second << "\t" << _curAmplicon->IsProperlyPaired();

    // update the report / log
    if (_report)
        _report << buffer.str() << std::endl;
    if (verbose)
        LOG_TRACE(buffer.str());
}

// _softmask performs the CIGAR string adjustment for the current record.
// if maskPrimers is true, primer sequence will also be softmasked.
void artic::Softmasker::_softmask(bool maskPrimers)
{

    // get the amplicon span, with or without primers
    std::pair<int64_t, int64_t> span = (maskPrimers) ? _curAmplicon->GetMinSpan() : _curAmplicon->GetMaxSpan();

    // update a counter before trimming
    if ((_curRec->core.pos < span.first) || (bam_endpos(_curRec) > span.second))
        _trimCounter++;

    // TODO: catch trim errors / report them via debug
    // if start of alignment is before amplicon start, mask
    if (_curRec->core.pos < span.first)
        TrimAlignment(_curRec, span.first, false);

    // if end of alignment is after amplicon end, mask
    if (bam_endpos(_curRec) > span.second)
        TrimAlignment(_curRec, span.second, true);
}

// Softmasker constructor.
artic::Softmasker::Softmasker(artic::PrimerScheme* primerScheme, const std::string& bamFile, const std::string& userCmd, unsigned int minMAPQ, unsigned int normalise, bool removeBadPairs, bool noReadGroups, bool primerStart, const std::string& reportFilename)
    : _primerScheme(primerScheme), _minMAPQ(minMAPQ), _normalise(normalise), _removeBadPairs(removeBadPairs), _noReadGroups(noReadGroups), _maskPrimerStart(primerStart)
{

    // get the input BAM or use STDIN if none given
    if (!bamFile.empty())
    {
        _inputBAM = sam_open(bamFile.c_str(), "r");
        if (!_inputBAM)
        {
            throw std::runtime_error("failed to open bam file: " + bamFile);
        }
    }
    else
    {
        _inputBAM = sam_open("-", "r");
        if (!_inputBAM)
        {
            throw std::runtime_error("cannot read BAM from STDIN - make sure you are piping a BAM file");
        }
    }

    // update the header with the called command and the primer pools
    _bamHeader = sam_hdr_read(_inputBAM);
    if (!_bamHeader)
        throw std::runtime_error("cannot access BAM header");
    artic::AddPGtoHeader(&_bamHeader, userCmd);
    if (!_noReadGroups)
    {
        for (auto pool : _primerScheme->GetPrimerPools())
            artic::AddRGtoHeader(&_bamHeader, pool);
    }

    // setup a report file if requested
    if (!reportFilename.empty())
    {
        _report.open(reportFilename, std::fstream::out | std::fstream::app);
        _report << "QueryName\tReferenceStart\tReferenceEnd\tPrimerPair\tPrimer1\tPrimer1Start\tPrimer2\tPrimer2Start\tIsSecondary\tIsSupplementary\tStart\tEnd\tCorrectlyPaired" << std::endl;
    }

    // get the holders ready
    _curRec = bam_init1();
    _recordCounter = 0;
    _filterDroppedCounter = 0;
    _normaliseDroppedCounter = 0;
    _trimCounter = 0;
}

// Softmasker destructor.
artic::Softmasker::~Softmasker(void)
{
    if (_curRec)
        bam_destroy1(_curRec);
    if (_bamHeader)
        bam_hdr_destroy(_bamHeader);
    if (_inputBAM)
        hts_close(_inputBAM);
    if (_report)
        _report.close();
}

// Run will perform the softmasking on the open BAM file.
void artic::Softmasker::Run(bool verbose)
{
    artic::Log::Init("softmasker");
    LOG_INFO("starting softmasking");
    if (_maskPrimerStart)
        LOG_INFO("include primers in softmask: true");

    // open up a new BAM for the output
    htsFile* outBam = hts_open("-", "wb");
    if (!outBam)
        throw std::runtime_error("cannot open BAM stream for writing");
    if (sam_hdr_write(outBam, _bamHeader) < 0)
        throw std::runtime_error("could not write header to output BAM stream");

    // iterate over the input BAM records
    while (sam_read1(_inputBAM, _bamHeader, _curRec) >= 0)
    {
        _recordCounter++;

        // skip unmapped and supplementary alignment records
        auto failedCheck = _checkRecord();
        if (failedCheck)
        {
            LOG_WARN("{} {}", bam_get_qname(_curRec), getErrorMsg(failedCheck));
            _filterDroppedCounter++;
            continue;
        }

        // get predicted amplicon for this alignment record based on the nearest primers
        auto amplicon = _primerScheme->FindPrimers(_curRec->core.pos, bam_endpos(_curRec));
        _curAmplicon = &amplicon;

        // add a primer pool readgroup to the alignment record based on the primer pairing
        // NOTE: primerscheme logic has already added "unmatched" as the primer pool if the FindPrimers method returns primers which are not properly paired
        if (!_noReadGroups)
            bam_aux_append(_curRec, "RG", 'Z', _curAmplicon->GetPrimerPool().size() + 1, (uint8_t*)_curAmplicon->GetPrimerPool().c_str());

        if (_removeBadPairs && !_curAmplicon->IsProperlyPaired())
        {
            LOG_WARN("{} skipped as not correctly paired ({})", bam_get_qname(_curRec), _curAmplicon->GetID());
            _filterDroppedCounter++;
            continue;
        }

        // if requested, update the report/stderr with this alignment record + amplicon details
        if (_report || verbose)
            _reportLine(verbose);

        //
        // TODO: collect k-mer frequencies
        //

        // stop processing the alignment record if normalise threshold reached for this amplicon
        if (_getAmpliconCount() >= _normalise)
        {
            LOG_WARN("{} dropped as abundance threshold reached", bam_get_qname(_curRec));
            _normaliseDroppedCounter++;
            continue;
        }

        // softmask amplicon, either to amplicon start or end
        _softmask(_maskPrimerStart);
        if (sam_write1(outBam, _bamHeader, _curRec) < 0)
            throw std::runtime_error("could not write record");
    }

    // print some stats
    LOG_INFO("finished softmasking");
    LOG_INFO("-\t{} alignments processed", _recordCounter);
    LOG_INFO("-\t{} alignments dropped by filters", _filterDroppedCounter);
    LOG_INFO("-\t{} alignments dropped after normalisation", _normaliseDroppedCounter);
    LOG_INFO("-\t{} alignments trimmed within amplicons", _trimCounter);
    if (verbose)
    {
        LOG_TRACE("amplicon\talignment count");
        for (auto const& amplicon : _ampliconCounter)
            LOG_TRACE("{}\t{}", amplicon.first, amplicon.second);
    }

    // close outfiles
    hts_close(outBam);
    return;
}

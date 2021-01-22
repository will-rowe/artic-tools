#include <algorithm>
#include <fstream>

#include "log.hpp"
#include "vcfCheck.hpp"
#include "version.hpp"

// VcfChecker constructor.
artic::VcfChecker::VcfChecker(artic::PrimerScheme* primerScheme, const std::string& vcfIn, const std::string& reportOut, const std::string& vcfOut, float minQual)
    : _primerScheme(primerScheme), _inputVCFfilename(vcfIn), _outputReportfilename(reportOut), _outputVCFfilename(vcfOut), _minQual(minQual)
{
    // get the input VCF ready
    _inputVCF = bcf_open(_inputVCFfilename.c_str(), "r");
    if (!_inputVCF)
        throw std::runtime_error("unable to open VCF file for reading");
    _vcfHeader = bcf_hdr_read(_inputVCF);
    if (!_vcfHeader)
        throw std::runtime_error("unable to read VCF header");
    _curRec = bcf_init();
    _prevRec = bcf_init();

    // get the outputs ready
    if (vcfOut.size() != 0)
    {
        _outputVCF = bcf_open(vcfOut.c_str(), "w");
        if (!_outputVCF)
            throw std::runtime_error("unable to open output VCF file for writing");
        std::string progLine = "##" + PROG_NAME + "_version=" + artic::GetVersion();
        bcf_hdr_append(_vcfHeader, progLine.c_str());
        if (bcf_hdr_write(_outputVCF, _vcfHeader) < 0)
            throw std::runtime_error("could not write header to output VCF stream");
    }
    else
    {
        _outputVCF = nullptr;
    }

    // zero counters
    _recordCounter = 0;
    _keepCounter = 0;
    _numValid = 0;
    _numPrimerSeq = 0;
    _numAmpOverlap = 0;
    _numAmpOverlapDiscard = 0;
    _numLowQual = 0;
}

// VcfChecker destructor.
artic::VcfChecker::~VcfChecker(void)
{
    if (_inputVCF)
        bcf_close(_inputVCF);
    if (_outputVCF)
        bcf_close(_outputVCF);
    if (_vcfHeader)
        bcf_hdr_destroy(_vcfHeader);
    if (_curRec)
        bcf_destroy(_curRec);
    if (_prevRec)
        bcf_destroy(_prevRec);
}

// checkRecordValidity will return true if the current record passes the basic validity checks.
bool artic::VcfChecker::_checkRecordValidity()
{
    // check var reference is in primer scheme
    std::string refID = bcf_hdr_id2name(_vcfHeader, _curRec->rid);
    if (refID != _primerScheme->GetReferenceName())
    {
        LOG_ERROR("\tdropping - reference ID does not match primer scheme reference ({})", refID);
        return false;
    }

    // check var primer pool is specified and is in scheme
    int ndst = 0;
    char* dst = NULL;
    bcf_get_info_string(_vcfHeader, _curRec, "Pool", &dst, &ndst);
    if (ndst == 0)
    {
        LOG_ERROR("\tdropping - no pool information provided");
        if (dst)
            free(dst);
        return false;
    }

    /*
        std::string pool(dst);
        std::vector<std::string>::iterator poolItr = std::find(primerPools.begin(), primerPools.end(), pool);
        if (poolItr == primerPools.end())
        {
            LOG_ERROR("\tdropping - pool not found in scheme ({})", pool);
            return false;
        }
    */

    if (dst)
        free(dst);

    // check var position is in scheme bounds
    if (_curRec->pos < _primerScheme->GetRefStart() || _curRec->pos > _primerScheme->GetRefEnd())
    {
        LOG_ERROR("\tdropping - outside of scheme bounds ({}:{})", _primerScheme->GetRefStart(), _primerScheme->GetRefEnd());
        return false;
    }
    return true;
}

//
void artic::VcfChecker::_getRecordStats()
{
    // keepOverlaps is used to write to record in one pass
    // if they are both identical vars in the same overlap pos.
    bool keepOverlaps = false;

    // evaluate if the prev record is to be kept or not
    bcf_unpack(_prevRec, BCF_UN_STR);
    auto adjustedPos = _prevRec->pos + 1;
    LOG_TRACE("variant at pos {}: {}->{}", adjustedPos, _prevRec->d.allele[0], _prevRec->d.allele[1]);
    bool discardRec = false;

    // check if in primer site
    if (_primerScheme->CheckPrimerSite(_prevRec->pos))
    {
        LOG_WARN("\tlocated within a primer sequence");
        _numPrimerSeq++;
    }

    // check quality
    if (_prevRec->qual < _minQual)
    {
        LOG_WARN("\tqual ({}) is below minimum quality threshold", _prevRec->qual);
        _numLowQual++;
        discardRec = true;
    }

    // check amplicon overlap
    if (_primerScheme->CheckAmpliconOverlap(_prevRec->pos))
    {
        LOG_TRACE("\tlocated within an amplicon overlap region");
        _numAmpOverlap++;
        if (_prevRec->pos != _curRec->pos)
        {
            // check the position of the next var in the file
            //LOG_ERROR("\tvar pos does not match with that of previously identified overlap var, holding new var (and dropping held var at {})", _prevRec->pos + 1);
            _numAmpOverlapDiscard++;
            discardRec = true;
        }

        else if (strcmp(_prevRec->d.allele[1], _curRec->d.allele[1]))
        {
            // check if the held var is the same allele
            // TODO: this logic won't check all against all if num allele > 2
            //LOG_ERROR("\tvar pos matches that of previously identified overlap var but alleles mismatch, holding new var (and dropping held var at {})", adjustedPos);
            _numAmpOverlapDiscard++;
            discardRec = true;
        }

        else
        {
            // held var matches the current var, so we can write the held record, clear the holder, and let the loop progress so that the current rec is also written
            // TODO: this would be the place to merge copies as discussed at https://github.com/will-rowe/artic-tools/issues/3
            //LOG_TRACE("\tmultiple copies of var found at pos {} in overlap region, keeping all copies", adjustedPos);
            keepOverlaps = true;
            _keepCounter++; // increment for one of the two records
        }
    }

    // if the record passes all checks, keep it
    if (!discardRec)
    {
        _keepCounter++;
        if (_outputVCF)
        {
            if (bcf_write(_outputVCF, _vcfHeader, _prevRec) < 0)
                throw std::runtime_error("could not write record");
            if (keepOverlaps)
            {
                if (bcf_write(_outputVCF, _vcfHeader, _curRec) < 0)
                    throw std::runtime_error("could not write record");
            }
        }
    }

    // update the previous record with the current record
    bcf_empty(_prevRec);
    _prevRec = bcf_dup(_curRec);
}

// Run will perform the softmasking on the open BAM file.
void artic::VcfChecker::Run()
{
    LOG_TRACE("setting parameters");
    LOG_TRACE("\toutput report: {}", _outputReportfilename);
    if (_outputVCF)
    {
        LOG_TRACE("\tfiltering variants: true");
        LOG_TRACE("\toutput VCF: {}", _outputVCFfilename);
    }
    else
    {
        LOG_TRACE("\tfiltering variants: false");
    }
    LOG_TRACE("\tminimum quality threshold: {}", _minQual);

    // iterate VCF file
    LOG_TRACE("reading VCF file");
    while (bcf_read(_inputVCF, _vcfHeader, _curRec) == 0)
    {

        // unpack the record and log it
        if (_curRec->errcode)
            throw std::runtime_error("refusing to process VCF records");
        bcf_unpack(_curRec, BCF_UN_STR);
        _recordCounter++;

        // check the current record is valid
        if (!_checkRecordValidity())
            continue;
        _numValid++;

        // if this is the first record, put it in the checking holder and move to the next record
        if (_numValid == 1)
        {
            _prevRec = bcf_dup(_curRec);
            continue;
        }

        // get the stats
        _getRecordStats();
    }

    // get stats from final record
    _getRecordStats();

    // write the stats to the report file
    std::ofstream outputReport;
    outputReport.open(_outputReportfilename);
    outputReport << "total vars.\tinvalid vars.\tin primer sites vars.\tin amplicon overlaps vars.\tqual <" << _minQual << " vars.\tinput VCF file" << std::endl;
    outputReport << _recordCounter << "\t" << _recordCounter - _numValid << "\t" << _numPrimerSeq << "\t" << _numAmpOverlap << "\t" << _numLowQual << "\t" << _inputVCFfilename << std::endl;
    outputReport.close();

    // finish up
    LOG_TRACE("finished checking")
    LOG_INFO("\t{} variant records processed", _recordCounter);
    LOG_INFO("\t{} variant records passed checks", _keepCounter);
}

// GetNumRecords returns the number of records processed by the checker.
unsigned int artic::VcfChecker::GetNumRecords(void) const { return _recordCounter; }

// GetNumInPrimerSite returns the number of records within primer sites.
unsigned int artic::VcfChecker::GetNumInPrimerSite(void) const { return _numPrimerSeq; }

// GetNumInOverlap returns the number of records within amplicon overlap regions.
unsigned int artic::VcfChecker::GetNumInOverlap(void) const { return _numAmpOverlap; }

// GetNumLowQual returns the number of records with low quality.
unsigned int artic::VcfChecker::GetNumLowQual(void) const { return _numLowQual; }
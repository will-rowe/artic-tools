#include <algorithm>
#include <fstream>

#include "log.hpp"
#include "vcfCheck.hpp"
#include "version.hpp"

// VcfChecker constructor.
artic::VcfChecker::VcfChecker(artic::PrimerScheme* primerScheme, const std::string& vcfIn, const std::string& reportOut, const std::string& vcfOut, bool dropOverlapFails, float minQual)
    : _primerScheme(primerScheme), _outputReportfilename(reportOut), _outputVCFfilename(vcfOut), _dropOverlapFails(dropOverlapFails), _minQual(minQual)
{
    // get the input VCF ready
    _inputVCF = bcf_open(vcfIn.c_str(), "r");
    if (!_inputVCF)
        throw std::runtime_error("unable to open VCF file for reading");
    _vcfHeader = bcf_hdr_read(_inputVCF);
    if (!_vcfHeader)
        throw std::runtime_error("unable to read VCF header");
    _curRec = bcf_init();
    _dupCheck = false;
    _recHolder = bcf_init();

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

    // zero counters
    _recordCounter = 0;
    _keepCounter = 0;
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
    if (_recHolder)
        bcf_destroy(_recHolder);
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

// Run will perform the softmasking on the open BAM file.
void artic::VcfChecker::Run()
{
    LOG_TRACE("setting parameters");
    LOG_TRACE("\toutput report: {}", _outputReportfilename);
    if (_outputVCF)
    {
        LOG_TRACE("\toutput VCF: {}", _outputVCFfilename);
        LOG_TRACE("\tfiltering variants: true");
        LOG_TRACE("\tdiscard overlap fail vars: {}", _dropOverlapFails);
    }
    else
    {
        LOG_TRACE("\tfiltering variants: false");
    }
    LOG_TRACE("\tminimum quality threshold: {}", _minQual);

    // vcf stats
    int numInvalid = 0;
    int numPrimerSeq = 0;
    int numAmpOverlap = 0;
    int numLowQual = 0;

    // iterate VCF file
    LOG_TRACE("reading VCF file");
    auto primerPools = _primerScheme->GetPrimerPools();
    std::string prevAl;
    bool qualDiscard = false;
    while (bcf_read(_inputVCF, _vcfHeader, _curRec) == 0)
    {

        // unpack the record and log it
        if (_curRec->errcode)
            throw std::runtime_error("refusing to process VCF records");
        bcf_unpack(_curRec, BCF_UN_STR);
        _recordCounter++;
        auto adjustedPos = _curRec->pos + 1;
        LOG_TRACE("variant at pos {}: {}->{}", adjustedPos, _curRec->d.allele[0], _curRec->d.allele[1]);

        // check the current record is valid
        if (!_checkRecordValidity())
        {
            numInvalid++;
            continue;
        }

        // check if in primer site
        if (_primerScheme->CheckPrimerSite(_curRec->pos))
        {
            LOG_WARN("\tlocated within a primer sequence");
            numPrimerSeq++;
        }

        // check quality
        if (_curRec->qual < _minQual)
        {
            LOG_WARN("\tqual ({}) is below minimum quality threshold", _curRec->qual);
            numLowQual++;
            continue;
        }

        // check amplicon overlap
        if (_primerScheme->CheckAmpliconOverlap(_curRec->pos))
        {
            LOG_TRACE("\tlocated within an amplicon overlap region");

            // check if we've already seen a var in this overlap position
            if (!_dupCheck)
            {
                LOG_TRACE("\tnothing seen at position yet, holding var");
                _recHolder = bcf_dup(_curRec);
                _dupCheck = true;
                prevAl = std::string(_curRec->d.allele[1]);
                continue;
            }

            // check the held var is at the same position
            if (_curRec->pos != _recHolder->pos)
            {
                LOG_ERROR("\tvar pos does not match with that of previously identified overlap var, holding new var (and dropping held var at {})", _recHolder->pos + 1);
                bcf_empty(_recHolder);
                _recHolder = bcf_dup(_curRec);
                prevAl = std::string(_curRec->d.allele[1]);
                numAmpOverlap++;
                if (!_dropOverlapFails && _outputVCF)
                {
                    _keepCounter++;
                    if (bcf_write(_outputVCF, _vcfHeader, _recHolder) < 0)
                        throw std::runtime_error("could not write record");
                }
                continue;
            }

            // now check if the held var is the same allele
            // TODO: this logic won't check all against all if num allele > 2
            if (strcmp(_curRec->d.allele[1], prevAl.c_str()))
            {
                LOG_ERROR("\tvar pos matches that of previously identified overlap var but alleles mismatch, holding new var (and dropping held var at {})", adjustedPos);
                bcf_empty(_recHolder);
                _recHolder = bcf_dup(_curRec);
                prevAl = std::string(_curRec->d.allele[1]);
                numAmpOverlap++;
                if (!_dropOverlapFails && _outputVCF)
                {
                    _keepCounter++;
                    if (bcf_write(_outputVCF, _vcfHeader, _recHolder) < 0)
                        throw std::runtime_error("could not write record");
                }
                continue;
            }

            // held var matches the current var, so we can write the held record, clear the holder, and let the loop progress so that the current rec is also written
            // TODO: this would be the place to merge copies as discussed at https://github.com/will-rowe/artic-tools/issues/3
            LOG_TRACE("\tmultiple copies of var found at pos {} in overlap region, keeping all copies", adjustedPos);
            if (_outputVCF)
                if (bcf_write(_outputVCF, _vcfHeader, _recHolder) < 0)
                    throw std::runtime_error("could not write record");
            _keepCounter++;
            bcf_empty(_recHolder);
            _dupCheck = false;
        }

        _keepCounter++;
        if (_outputVCF)
            if (bcf_write(_outputVCF, _vcfHeader, _curRec) < 0)
                throw std::runtime_error("could not write record");
    }

    // process anything left in the overlap check holder
    if (_dupCheck)
    {
        LOG_ERROR("\tdropping var at pos {} which is in an amplicon overlap region but was only found once", _recHolder->pos + 1);
        numAmpOverlap++;
        if (!_dropOverlapFails && _outputVCF)
        {
            _keepCounter++;
            if (bcf_write(_outputVCF, _vcfHeader, _recHolder) < 0)
                throw std::runtime_error("could not write record");
        }
    }

    // write the stats to the report file
    std::ofstream outputReport;
    outputReport.open(_outputReportfilename);
    outputReport << "num. invalid vars.\tvars. in primer sites\tvars. once in amplicon overlaps\tvars. qual <" << _minQual << std::endl;
    outputReport << numInvalid << "\t" << numPrimerSeq << "\t" << numAmpOverlap << "\t" << numLowQual << std::endl;
    outputReport.close();

    // finish up
    LOG_TRACE("finished checking")
    LOG_INFO("\t{} variant records processed", _recordCounter);
    LOG_INFO("\t{} variant records passed checks", _keepCounter);
}

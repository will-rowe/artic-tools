#include <algorithm>
#include <iostream>

#include "log.hpp"
#include "vcfCheck.hpp"
#include "version.hpp"

// VcfChecker constructor.
artic::VcfChecker::VcfChecker(artic::PrimerScheme* primerScheme, const std::string& vcfIn, const std::string& vcfOut, bool dropPrimerVars)
    : _primerScheme(primerScheme), _outfileName(vcfOut), _dropPrimerVars(dropPrimerVars)
{
    // get the input VCF ready
    _inputVCF = bcf_open(vcfIn.c_str(), "r");
    if (!_inputVCF)
        throw std::runtime_error("unable to open VCF file for reading");
    _vcfHeader = bcf_hdr_read(_inputVCF);
    if (!_vcfHeader)
        throw std::runtime_error("unable to read VCF header");
    _curRec = bcf_init();

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
}

// Run will perform the softmasking on the open BAM file.
void artic::VcfChecker::Run(bool noLog)
{
    if (!noLog)
    {
        artic::Log::Init("vcfchecker");
        LOG_INFO("starting VCF checker");
    }

    // get the output file ready if neeeded
    if ((_outfileName.size() != 0))
    {
        _outputVCF = bcf_open(_outfileName.c_str(), "w");
        if (!_outputVCF)
            throw std::runtime_error("unable to open VCF file for writing");
        if (!noLog)
            LOG_INFO("writing variants passing checks to {}", _outfileName);

        // add prog info and copy the header
        std::string progLine = "##" + PROG_NAME + "_version=" + artic::GetVersion();
        bcf_hdr_append(_vcfHeader, progLine.c_str());
        if (bcf_hdr_write(_outputVCF, _vcfHeader) < 0)
            throw std::runtime_error("could not write header to output VCF stream");
    }

    // holders
    auto primerPools = _primerScheme->GetPrimerPools();
    int ndst = 0;
    char* dst = NULL;

    // iterate VCF file
    while (bcf_read(_inputVCF, _vcfHeader, _curRec) == 0)
    {
        _recordCounter++;

        // check var reference is in primer scheme
        std::string refID = bcf_hdr_id2name(_vcfHeader, _curRec->rid);
        if (refID != _primerScheme->GetReferenceID())
        {
            if (!noLog)
                LOG_WARN("variant reference ID does not match primer scheme reference - {}", refID);
            continue;
        }

        // check var position is in scheme bounds
        if (_curRec->pos < _primerScheme->GetRefStart() || _curRec->pos > _primerScheme->GetRefEnd())
        {
            if (!noLog)
                LOG_WARN("variant outside of scheme bounds - {} at {}", _curRec->d.allele[1], _curRec->pos);
            continue;
        }

        // check var primer pool is specified and is in scheme
        bcf_get_info_string(_vcfHeader, _curRec, "Pool", &dst, &ndst);
        if (ndst == 0)
        {
            if (!noLog)
                LOG_WARN("no Pool information for variant at {}", _curRec->pos);
            continue;
        }
        std::string pool(dst);
        std::vector<std::string>::iterator poolItr = std::find(primerPools.begin(), primerPools.end(), pool);
        if (poolItr == primerPools.end())
        {
            if (!noLog)
                LOG_WARN("variant Pool not in scheme - {}", dst);
            continue;
        }

        // check if in primer site
        if (_primerScheme->CheckPrimerSite(_curRec->pos, pool))
        {
            if (!noLog)
                LOG_WARN("variant found within primer sequence - {}", _curRec->pos);
            if (_dropPrimerVars)
                continue;
        }

        // check amplicon overlap
        if (_primerScheme->CheckAmpliconOverlap(_curRec->pos))
        {
            if (!noLog)
                LOG_WARN("variant found in scheme region with amplicon overlap - {}", _curRec->pos);

            // add to a list
        }

        _keepCounter++;

        // TODO: optional - add in depth, qual filters etc.
        //bcf_unpack(_curRec, BCF_UN_STR);
        //std::cout << _curRec->pos << " " << _curRec->rlen << " " << seqnames[_curRec->rid] << std::endl;
        //LOG_INFO("variant okay: {} -> {} at {}", _curRec->d.allele[0], _curRec->d.allele[1], _curRec->pos);
    }

    // free holders
    if (dst)
        free(dst);

    // print some stats
    if (!noLog)
    {
        LOG_INFO("{} variant records processed", _recordCounter);
        LOG_INFO("{} variant records passed checks", _keepCounter);
    }
}

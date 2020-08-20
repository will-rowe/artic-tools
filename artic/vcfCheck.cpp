#include <algorithm>
#include <iostream>

#include "log.hpp"
#include "vcfCheck.hpp"
#include "version.hpp"

// VcfChecker constructor.
artic::VcfChecker::VcfChecker(artic::PrimerScheme* primerScheme, const std::string& vcfIn, const std::string& vcfOut, bool dropPrimerVars, bool dropOverlapFails)
    : _primerScheme(primerScheme), _outfileName(vcfOut), _dropPrimerVars(dropPrimerVars), _dropOverlapFails(dropOverlapFails)
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

// Run will perform the softmasking on the open BAM file.
void artic::VcfChecker::Run()
{
    std::cout << _outfileName << std::endl;
    artic::Log::Init("vcfchecker");
    LOG_TRACE("starting VCF checker");
    if ((_outfileName.size() != 0))
    {
        LOG_TRACE("\tfiltering variants: true");
        LOG_TRACE("\toutput file: {}", _outfileName);
    }
    else
    {
        _outputVCF = 0;
        LOG_TRACE("\tfiltering variants: false");
    }
    LOG_TRACE("\tdiscard primer site vars: {}", _dropPrimerVars);
    LOG_TRACE("\tdiscard overlap fail vars: {}", _dropOverlapFails);

    // get the output file ready if neeeded
    if ((_outfileName.size() != 0))
    {
        _outputVCF = bcf_open(_outfileName.c_str(), "w");
        if (!_outputVCF)
            throw std::runtime_error("unable to open VCF file for writing");
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
        bcf_unpack(_curRec, BCF_UN_STR);
        _recordCounter++;
        auto adjustedPos = _curRec->pos + 1;
        LOG_TRACE("variant at pos {}: {}->{}", adjustedPos, _curRec->d.allele[0], _curRec->d.allele[1]);

        // check var reference is in primer scheme
        std::string refID = bcf_hdr_id2name(_vcfHeader, _curRec->rid);
        if (refID != _primerScheme->GetReferenceName())
        {
            LOG_ERROR("\tdropping - reference ID does not match primer scheme reference ({})", refID);
            continue;
        }

        // check var primer pool is specified and is in scheme
        bcf_get_info_string(_vcfHeader, _curRec, "Pool", &dst, &ndst);
        if (ndst == 0)
        {
            LOG_ERROR("\tdropping - no pool information provided");
            continue;
        }
        std::string pool(dst);
        std::vector<std::string>::iterator poolItr = std::find(primerPools.begin(), primerPools.end(), pool);
        if (poolItr == primerPools.end())
        {
            LOG_ERROR("\tdropping - pool not found in scheme ({})", pool);
            continue;
        }

        // check var position is in scheme bounds
        if (_curRec->pos < _primerScheme->GetRefStart() || _curRec->pos > _primerScheme->GetRefEnd())
        {
            LOG_ERROR("\tdropping - outside of scheme bounds ({}:{})", _primerScheme->GetRefStart(), _primerScheme->GetRefEnd());
            continue;
        }

        // check if in primer site
        if (_primerScheme->CheckPrimerSite(_curRec->pos, pool))
        {
            if (_dropPrimerVars)
            {
                LOG_ERROR("\tdropping - located within a primer sequence for the primer pool ({})", pool);
                continue;
            }
            LOG_WARN("\tlocated within a primer sequence for the primer pool ({})", pool);
        }

        // check amplicon overlap
        // todo: handle if multiple vars found at a pos but alleles are different?
        // todo: merge identical vars in overlaps?
        if (_primerScheme->CheckAmpliconOverlap(_curRec->pos))
        {
            LOG_TRACE("\tlocated within an amplicon overlap region");

            // check if we've already seen an var in this overlap position
            if (!_dupCheck)
            {
                LOG_TRACE("\tnothing seen at position yet, holding var");
                _recHolder = bcf_dup(_curRec);
                _dupCheck = true;
                continue;
            }
            if (_curRec->pos != _recHolder->pos)
            {
                LOG_ERROR("\tvar pos does not match with that of previously identified overlap, holding var (and dropping held var at {})", _recHolder->pos);
                bcf_empty(_recHolder);
                _recHolder = bcf_dup(_curRec);
                continue;
            }

            // otherwise, the write the record we had on hold and clear the holder
            // TODO: this would be the place to merge copies as discussed at https://github.com/will-rowe/artic-tools/issues/3
            LOG_TRACE("\tmultiple copies of var found at pos {} in overlap region, keeping all copies", adjustedPos);
            if (_outfileName.size() != 0)
                if (bcf_write(_outputVCF, _vcfHeader, _recHolder) < 0)
                    throw std::runtime_error("could not write record");
            _keepCounter++;
            bcf_empty(_recHolder);
            _dupCheck = false;
        }
        _keepCounter++;
        if (_outfileName.size() != 0)
            if (bcf_write(_outputVCF, _vcfHeader, _curRec) < 0)
                throw std::runtime_error("could not write record");
    }

    // free holders
    if (dst)
        free(dst);

    // finish up
    LOG_TRACE("finished checking")
    if (_dupCheck)
    {
        LOG_ERROR("\tdropped var at {} which is in an amplicon overlap region but was only found once", _recHolder->pos);
    }
    LOG_INFO("\t{} variant records processed", _recordCounter);
    LOG_INFO("\t{} variant records passed checks", _keepCounter);
}

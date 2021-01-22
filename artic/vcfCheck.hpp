#ifndef VCF_FILTER_H
#define VCF_FILTER_H

#include <htslib/vcf.h>

#include "primerScheme.hpp"

namespace artic
{

    //******************************************************************************
    // VcfChecker class handles the VCf filtering
    //******************************************************************************
    class VcfChecker
    {
    public:
        // VcfChecker constructor and destructor.
        VcfChecker(artic::PrimerScheme* primerScheme, const std::string& vcfIn, const std::string& reportOut, const std::string& vcfOut, float minQual);
        ~VcfChecker(void);

        // Run will perform the filtering on the open VCF file.
        void Run();

    private:
        bool _checkRecordValidity(); // checkRecordValidity will return true if the current record passes the basic validity checks
        void _getRecordStats();

        // data holders
        artic::PrimerScheme* _primerScheme; // the loaded primer scheme
        std::string _outputReportfilename;  // the filename of the output report
        std::string _outputVCFfilename;     // the filename of the output VCF file
        vcfFile* _inputVCF;                 // the input VCF for filtering
        bcf_hdr_t* _vcfHeader;              // the input VCF header
        bcf1_t* _curRec;                    // the current VCF record being processed
        bcf1_t* _prevRec;                   // the previous VCF record, used to check against
        vcfFile* _outputVCF;                // the output VCF file

        // user parameters
        float _minQual; // the QUAL threshold for keeping records

        // counters
        unsigned int _recordCounter; // number of records processed by the softmasker
        unsigned int _keepCounter;   // number of records which passed filters
        unsigned int _numValid;
        unsigned int _numPrimerSeq;
        unsigned int _numAmpOverlap;
        unsigned int _numAmpOverlapDiscard;
        unsigned int _numLowQual;
    };

} // namespace artic

#endif
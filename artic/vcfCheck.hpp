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
        void Run(void);
        unsigned int GetNumRecords(void) const;      // returns the number of records processed by the checker
        unsigned int GetNumInPrimerSite(void) const; // returns the number of records within primer sites
        unsigned int GetNumInOverlap(void) const;    // returns the number of records within amplicon overlap regions
        unsigned int GetNumLowQual(void) const;      // returns the number of records with low quality

    private:
        bool _checkRecordValidity(); // checkRecordValidity will return true if the current record passes the basic validity checks
        void _getRecordStats();

        // data holders
        artic::PrimerScheme* _primerScheme; // the loaded primer scheme
        std::string _inputVCFfilename;      // the filename of the input VCF file
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
        unsigned int _recordCounter;        // number of records processed by the softmasker
        unsigned int _keepCounter;          // number of records which passed filters
        unsigned int _numValid;             // number of records which are considered valid (have primer pool info, in scheme boundary etc.)
        unsigned int _numPrimerSeq;         // number of records within a primer sequence
        unsigned int _numAmpOverlap;        // number of records within an amplicon overlap region (excludes primer seq)
        unsigned int _numAmpOverlapDiscard; // _numAmpOverlap - those which don't have 2 copies per site
        unsigned int _numLowQual;           // number of records which are below the set quality threshold
    };

} // namespace artic

#endif
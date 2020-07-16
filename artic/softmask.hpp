#ifndef SOFTMASK_H
#define SOFTMASK_H

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <string>

#include <artic/primerScheme.hpp>

namespace artic
{

    //******************************************************************************
    // Softmasker class handles the alignment softmasking
    //******************************************************************************
    class Softmasker
    {
    public:
        // Softmasker constructor and destructor.
        Softmasker(artic::PrimerScheme* primerScheme, const std::string& bamFile, const std::string& userCmd, unsigned int minMAPQ, unsigned int normalise, bool removeBadPairs, bool noReadGroups, const std::string& reportFilename);
        ~Softmasker(void);

        // Run will perform the softmasking on the open BAM file.
        void Run(const std::string& outFilePrefix, bool verbose);

    private:
        unsigned int _checkRecord(void);      // returns an error if the currently held record fails filters and should be skipped
        unsigned int _getAmpliconCount(void); // returns the number of times the queried amplicon has been seen before
        void _reportLine(bool verbose);       // add the primer information for the current record to the open report
        void _softmask(bool maskPrimers);     // performs the CIGAR string adjustment for the current record

        // data holders
        artic::PrimerScheme* _primerScheme; // the loaded primer scheme
        htsFile* _inputBAM;                 // the input BAM for softmasking
        bam_hdr_t* _bamHeader;              // the input BAM header
        bam1_t* _curRec;                    // the current alignment record being processed
        artic::Amplicon* _curAmplicon;      // the amplicon for the current alignment record
        std::fstream _report;               // the report file

        // user parameters
        unsigned int _minMAPQ;   // the MAPQ threshold for keeping records
        unsigned int _normalise; // the normalise threshold (set to 0 if normalisation not required)
        bool _removeBadPairs;    // ignore records where primers are incorrectly paired
        bool _noReadGroups;      // don't use read group information during soft masking

        // counters
        std::unordered_map<std::string, unsigned int> _ampliconCounter; // counts the amplicon pairs encountered during softmasking
        unsigned int _recordCounter;                                    // number of records processed by the softmasker
        unsigned int _filterDroppedCounter;                             // number of records which failed filters
        unsigned int _normaliseDroppedCounter;                          // number of records that were dropped post normalisation
        unsigned int _trimCounter;                                      // number of records that were trimmed within amplicon (either forward or reverse)
        unsigned int _ptrimCounter;                                     // number of records that were trimmed within amplicon primer sequences (either forward or reverse)
    };

}; // namespace artic

#endif
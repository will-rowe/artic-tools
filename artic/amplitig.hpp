#ifndef AMPLITIG_H
#define AMPLITIG_H

#include <mutex>
#include <shared_mutex>

#include "kmers.hpp"
#include "primerScheme.hpp"

namespace artic
{

    //******************************************************************************
    // Amplitigger class handles the amplicon read binning
    //******************************************************************************
    class Amplitigger
    {
    public:
        // Amplitigger constructor and destructor.
        Amplitigger(artic::PrimerScheme* primerScheme, const std::string& refFile, const std::vector<std::string> inputFiles, const std::string& userCmd, unsigned int kmerSize);
        //~Amplitigger(void);

        // Run will perform the softmasking on the open BAM file.
        void Run(bool verbose);

    private:
        // data holders
        artic::PrimerScheme* _primerScheme;                                          // the loaded primer scheme
        const std::string _refFile;                                                  // the reference fasta file
        const std::vector<std::string> _inputFiles;                                  // input FASTQ files
        std::unordered_map<artic::kmer_t, std::vector<unsigned int>> _primerKmerMap; // map of primer k-mers and their origins
        // std::unordered_map<std::string, artic::Amplicon> _amplicons; // map of amplicons

        mutable std::shared_mutex _mutex;

        // user parameters
        unsigned int _kmerSize;       // the k-mer size to use
        unsigned int _minPrimerKmers; // the minimum proportion of matching primer k-mers required for amplicon assignment

        // counters
        unsigned int _recordCounter;        // number of records processed by the softmasker
        unsigned int _filterDroppedCounter; // number of records which failed filters
    };

} // namespace artic

#endif
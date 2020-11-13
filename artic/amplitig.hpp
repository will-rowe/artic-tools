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
        Amplitigger(artic::PrimerScheme* primerScheme, const std::string& refFile, const std::vector<std::string> inputFiles, unsigned int kmerSize);
        //~Amplitigger(void);

        // Run will perform the amplitigging.
        void Run();

    private:
        // data holders
        artic::PrimerScheme* _primerScheme;         // the loaded primer scheme
        const std::string _refFile;                 // the reference fasta file
        const std::vector<std::string> _inputFiles; // input FASTQ files
        kmermap_t _primerKmerMap{};                 // map of primer k-mers and their origins
        // std::unordered_map<std::string, artic::Amplicon> _amplicons; // map of amplicons

        mutable std::shared_mutex _mutex;

        // user parameters
        // TODO: implement these as options
        unsigned int _kmerSize; // the k-mer size to use
        float _minPrimerKmers;  // the minimum proportion of matching primer k-mers required for amplicon assignment
        int _minReadLength;     // drop reads shorter than this length
        int _maxReadLength;     // drop reads longer than this length (default is to use max amplicon span in the scheme + 10%)

        // counters
        unsigned int _readCounter;  // number of reads processed by the softmasker
        unsigned int _droppedLong;  // number of reads which were deemed too long
        unsigned int _droppedShort; // number of reads which were demmed too short
    };

} // namespace artic

#endif
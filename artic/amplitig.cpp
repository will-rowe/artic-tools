#include "amplitig.hpp"
#include "fastqParser.hpp"
#include "kmers.hpp"
#include "log.hpp"

// Amplitigger constructor.
artic::Amplitigger::Amplitigger(artic::PrimerScheme* primerScheme, const std::string& refFile, const std::vector<std::string> inputFiles, unsigned int kmerSize)
    : _primerScheme(primerScheme), _refFile(refFile), _inputFiles(inputFiles), _kmerSize(kmerSize)
{

    // check the params
    if (_kmerSize > MAX_K_SIZE)
        throw std::runtime_error("requested k-mer size greater than maximum allowed size (" + std::to_string(MAX_K_SIZE) + ")");
    if (_kmerSize > _primerScheme->GetMinPrimerLen())
        throw std::runtime_error("requested k-mer size greater than the smallest primer in scheme (" + std::to_string(_primerScheme->GetMinPrimerLen()) + ")");
    if (_inputFiles.size() == 0)
        throw std::runtime_error("no FASTQ files provided");

    // get the holders ready
    _recordCounter = 0;
    _filterDroppedCounter = 0;

    // get the containers for the expected amplicons
    //for (auto amplicon : _primerScheme->GetExpAmplicons())
    //    _amplicons.emplace(amplicon.GetName(), amplicon);

    // load the primer k-mers
    LOG_TRACE("collecting primer k-mers");
    LOG_TRACE("\tk-mer size used:\t{}", _kmerSize);
    LOG_TRACE("\treference fasta file:\t{}", _refFile);
    _primerScheme->GetPrimerKmers(_refFile, _kmerSize, _primerKmerMap);
    LOG_TRACE("\ttotal unique k-mers:\t{}", _primerKmerMap.size());
}

// Run will perform the amplicon read binning on the open FASTQ file.
void artic::Amplitigger::Run(bool verbose)
{

    // get the FASTQ parser ready
    size_t fileNum = 0;
    std::string seq;
    artic::FastqReader fastqReader(_inputFiles);

    // get a k-mer holder
    artic::kmerset_t kmers;

    // get the read k-mers
    LOG_TRACE("collecting read k-mers");
    while (fastqReader.GetRecord(seq, fileNum) >= 0)
    {
        _recordCounter++;
        artic::GetEncodedKmers(seq.c_str(), seq.size(), _kmerSize, kmers);
        LOG_TRACE("\tgot {} kmers", kmers.size());
        kmers.clear();
    }
    fastqReader.Close();

    // print some stats
    LOG_TRACE("\tread processed:\t{}", _recordCounter);
}

/*
grab a read
get vector of k-mers
loop through the vector and check against primer k-mers
add vector of k-mers to k-mer freq map of matched amplicon

requires 3x loops:
once to get k-mers from sequence
once to check k-mers against primer k-mers
once to add k-mers to correct amplicon map

improvements for effiency:
heuristic to skip regions during primer check (i.e, found one primer so can skip ~100 k-mers to next likley place a primer is found)
check k-mer against map as it is generated - might be a better way

additions:
bloom filter to each amplicon to remove unique k-mers
remove/exclude primer k-mers (but how to clip reads before/after primer sequence to within amplicons?)

*/
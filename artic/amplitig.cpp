#include "amplitig.hpp"
#include "fastqParser.hpp"
#include "kmers.hpp"
#include "log.hpp"

// Amplitigger constructor.
artic::Amplitigger::Amplitigger(artic::PrimerScheme* primerScheme, const std::string& refFile, const std::vector<std::string> inputFiles, const std::string& userCmd, unsigned int kmerSize)
    : _primerScheme(primerScheme), _refFile(refFile), _inputFiles(inputFiles), _kmerSize(kmerSize)
{

    // check the params
    if (_kmerSize > MAX_K_SIZE)
        throw std::runtime_error("requested k-mer size greater than maximum allowed size (" + std::to_string(MAX_K_SIZE) + ")");
    if (_kmerSize > _primerScheme->GetMinPrimerLen())
        throw std::runtime_error("requested k-mer size greater than the smallest primer in scheme (" + std::to_string(_primerScheme->GetMinPrimerLen()) + ")");

    // get the holders ready
    //_curRec = bam_init1();
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
    LOG_INFO("blam");

    size_t fileNum = 0;
    std::string seq;
    FastqFile fastqReader(_inputFiles);

    while (fastqReader.GetRecord(seq, fileNum) >= 0)
    {
        LOG_INFO(seq);
    }
    fastqReader.Close();
}
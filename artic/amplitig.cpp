#include <algorithm>
#include <experimental/filesystem>
#include <iostream>
#include <kseq++/seqio.hpp>
#include <utility>

#include "amplitig.hpp"
#include "kmers.hpp"
#include "log.hpp"

using namespace klibpp;

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

    // get thresholds
    _minPrimerKmers = 0.9;
    _minReadLength = 100;
    _maxReadLength = _primerScheme->GetMaxAmpliconSpan() + (0.10 * _primerScheme->GetMaxAmpliconSpan());

    // get the holders ready
    _readCounter = 0;
    _droppedLong = 0;
    _droppedShort = 0;

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
void artic::Amplitigger::Run()
{
    // get some holders ready
    artic::kmerset_t kmers;
    std::vector<unsigned int> ampliconIDs;
    KSeq record;
    int seqLen;

    LOG_TRACE("processing");
    for (auto file : _inputFiles)
    {
        if (!std::filesystem::exists(file))
            throw std::runtime_error("supplied file does not exist:\t" + file);
        LOG_TRACE("\treading file:\t{}", file);

        SeqStreamIn iss(file.c_str());
        while (iss >> record)
        {
            seqLen = record.seq.size();

            // check read length
            _readCounter++;
            if (seqLen < _minReadLength)
            {
                _droppedShort++;
                continue;
            }
            if (seqLen > _maxReadLength)
            {
                _droppedLong++;
                continue;
            }

            // get the read k-mers
            artic::GetEncodedKmers(record.seq.c_str(), seqLen, _kmerSize, kmers);

            // assign read to amplicon
            for (auto kmer : kmers)
            {

                // check if read k-mer is matched to a primer k-mer
                auto it = _primerKmerMap.find(kmer);
                if (it != _primerKmerMap.end())
                {
                    ampliconIDs.reserve(ampliconIDs.size() + it->second.size());
                    ampliconIDs.insert(ampliconIDs.end(), it->second.begin(), it->second.end());
                }
            }
            sort(ampliconIDs.begin(), ampliconIDs.end());

            std::vector<std::pair<unsigned int, int>> ampliconCandidates;
            int chain = 0;
            for (size_t i = 1; i < ampliconIDs.size(); i++)
            {
                // if successive match, keep building the chain
                if (ampliconIDs.at(i - 1) == ampliconIDs.at(i))
                {
                    chain++;
                    continue;
                }

                // ignore empty chains
                if (chain == 0)
                    continue;

                // if first chain, add it to candidates and continue
                if (ampliconCandidates.empty())
                {
                    ampliconCandidates.emplace_back(std::make_pair(ampliconIDs.at(i - 1), chain));
                    chain = 0;
                    continue;
                }

                // if new chain is < than existing max, continue
                if (chain < ampliconCandidates.back().second)
                {
                    chain = 0;
                    continue;
                }

                // if new chain > existing max, replace and continue
                if (chain > ampliconCandidates.back().second)
                {
                    ampliconCandidates.pop_back();
                    ampliconCandidates.emplace_back(std::make_pair(ampliconIDs.at(i - 1), chain));
                    chain = 0;
                    continue;
                }

                // final opt, same value as existing max so add it as well
                ampliconCandidates.emplace_back(std::make_pair(ampliconIDs.at(i - 1), chain));
                chain = 0;
            }

            // process the likely amplicons
            switch (ampliconCandidates.size())
            {
                case 0:
                    //LOG_ERROR("no amplicon found....");
                    break;
                case 1:
                    LOG_TRACE("{}\t{}\t{}", record.name, _primerScheme->GetAmpliconName(ampliconCandidates.front().first), ampliconCandidates.front().second);

                    // compare biggest chain to max %k-mer content
                    // remove chains with multiple identities
                    break;
                default:
                    //LOG_ERROR("too many amplicons found - {} amplicons with {} k-mers", ampliconCandidates.size(), ampliconCandidates.front().second);
                    break;
            }

            // clear the holders ready for the next read
            kmers.clear();
            ampliconIDs.clear();
        }
    }

    // print some stats
    LOG_TRACE("finished processing reads")
    LOG_TRACE("\ttotal input reads:\t{}", _readCounter);
    LOG_TRACE("\ttotal dropped reads:\t{}", _droppedLong + _droppedShort);
    LOG_TRACE("\t- short reads (<{}):\t{}", _minReadLength, _droppedShort);
    LOG_TRACE("\t- long reads (>{}):\t{}", _maxReadLength, _droppedLong);
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
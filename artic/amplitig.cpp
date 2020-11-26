#include <algorithm>
#include <boost/filesystem.hpp>
#include <iostream>
#include <kseq++/seqio.hpp>
#include <utility>

#include "amplitig.hpp"
#include "kmers.hpp"
#include "log.hpp"

using namespace klibpp;

// Amplitigger constructor.
artic::Amplitigger::Amplitigger(artic::PrimerScheme* primerScheme, const std::string& refFile, const std::vector<std::string> inputFiles, unsigned int kmerSize, float kmerMatch)
    : _primerScheme(primerScheme), _refFile(refFile), _inputFiles(inputFiles), _kmerSize(kmerSize), _minPrimerKmers(kmerMatch)
{

    // check the params
    if (_kmerSize > MAX_K_SIZE)
        throw std::runtime_error("requested k-mer size greater than maximum allowed size (" + std::to_string(MAX_K_SIZE) + ")");
    if (_kmerSize > _primerScheme->GetMinPrimerLen())
        throw std::runtime_error("requested k-mer size greater than the smallest primer in scheme (" + std::to_string(_primerScheme->GetMinPrimerLen()) + ")");
    if (_inputFiles.size() == 0)
        throw std::runtime_error("no FASTQ files provided");

    // TODO: collect thresholds from user?
    _minReadLength = 100;
    _maxReadLength = _primerScheme->GetMaxAmpliconSpan() + (0.10 * _primerScheme->GetMaxAmpliconSpan());

    // get the holders ready
    _readCounter = 0;
    _droppedLong = 0;
    _droppedShort = 0;
    _droppedUnbinned = 0;
    _multibinned = 0;

    // load the primer k-mers
    LOG_TRACE("collecting primer k-mers");
    LOG_TRACE("\tk-mer size used:\t{}", _kmerSize);
    LOG_TRACE("\tk-mer matches required:\t{}%", _minPrimerKmers);
    LOG_TRACE("\treference fasta file:\t{}", _refFile);
    _primerScheme->GetPrimerKmers(_refFile, _kmerSize, _primerKmerMap);
    LOG_TRACE("\ttotal distinct k-mers:\t{}", _primerKmerMap.size());

    // TODO: filter k-mers and retain only mutually exclusive k-mers
    int meCount = 0;
    for (auto i : _primerKmerMap)
    {
        if (i.second.size() > 1)
        {
            _primerKmerMap.erase(i.first);
        }
        else
        {
            meCount++;
        }
    }
    LOG_TRACE("\ttotal mutually exclusive k-mers:\t{}", meCount);
    if (meCount != _primerKmerMap.size())
        throw std::runtime_error("did not remove k-mers...");
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
        if (!boost::filesystem::exists(file))
            throw std::runtime_error("supplied file does not exist:\t" + file);
        LOG_TRACE("\treading file:\t{}", file);

        SeqStreamIn iss(file.c_str());
        while (iss >> record)
        {
            _readCounter++;

            // check read length
            seqLen = record.seq.size();
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

            // clear sets
            kmers.clear();
            ampliconIDs.clear();

            // get the read k-mers
            artic::GetEncodedKmers(record.seq.c_str(), seqLen, _kmerSize, kmers);

            // check read k-mers against primer scheme k-mers, keep amplicon IDs for all matches
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

            // sort the matching IDs and then find the best candidate amplicon ID for this read
            std::sort(ampliconIDs.begin(), ampliconIDs.end());
            std::vector<std::pair<unsigned int, int>> ampliconCandidates;
            int chain = 0;
            for (size_t i = 1; i < ampliconIDs.size(); i++)
            {
                // if successive match, continue to the next ID and keep building the chain
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
            int binned = 0;
            for (auto candidate : ampliconCandidates)
            {
                auto amplicon = _primerScheme->GetAmplicon(candidate.first);
                auto ampliconKmers = (amplicon.GetForwardPrimer()->GetLen() + amplicon.GetReversePrimer()->GetLen()) - (_kmerSize * 2) + 2;
                auto propKmers = float(candidate.second) / float(ampliconKmers);
                if (propKmers >= _minPrimerKmers)
                {
                    std::cout << record.name << "\t" << amplicon.GetName() << "\t" << propKmers << std::endl;
                    binned++;
                }
            }

            // update some numbers
            if (binned > 1)
                _multibinned++;

            // this will catch reads which didn't get any primer k-mer hits AND those which had too small a match proportion
            if (binned == 0)
                _droppedUnbinned++;

            // compare biggest chain to max %k-mer content
            // remove chains with multiple identities
        }
    }

    // print some stats
    LOG_TRACE("finished processing reads")
    LOG_TRACE("\ttotal input reads:\t{}", _readCounter);
    LOG_TRACE("\ttotal dropped reads:\t{}", (_droppedLong + _droppedShort + _droppedUnbinned));
    LOG_TRACE("\t- short reads (<{}):\t{}", _minReadLength, _droppedShort);
    LOG_TRACE("\t- long reads (>{}):\t{}", _maxReadLength, _droppedLong);
    LOG_TRACE("\t- unbinned reads:\t{}", _droppedUnbinned);
    LOG_TRACE("\ttotal binned reads:\t{}", (_readCounter - (_droppedLong + _droppedShort + _droppedUnbinned)));
    LOG_TRACE("\t- multibinned reads:\t{}", _multibinned);
}

/*
collect primer k-mers
remove any k-mers shared between primers
check each primer has some primers left / select a minimizer?

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
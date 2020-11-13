#include <iostream>
#include <queue>

#include "kmers.hpp"

namespace artic
{

    // GetEncodedKmers will compute and integer encode all canonical k-mers in a sequence, adding them to the provided container.
    void GetEncodedKmers(const char* seq, uint32_t seqLen, uint32_t kSize, kmerset_t& kmers)
    {
        if (kSize > MAX_K_SIZE)
            throw std::runtime_error("k-mer size must be <= " + std::to_string(MAX_K_SIZE));
        //if (!kmers.empty())
        //    throw std::runtime_error("provided k-mer container already has k-mers in it");
        uint64_t kmerMask = (1ULL << (2 * kSize)) - 1;
        uint64_t bitShift = 2 * (kSize - 1);
        kmer_t x[2];
        x[0] = x[1] = 0;
        uint32_t l = 0;
        for (uint32_t i = 0; i < seqLen; i++)
        {

            // get 2-bit encoding of current base
            uint8_t base = (uint8_t)seq[i] < 128 ? nt2char[static_cast<uint8_t>(seq[i])] : 4;

            // skip non A/C/T/G bases
            if (base > 3)
                continue;

            // use the masks to add the current base to the k-mer and its reverse complement, dropping the oldest base from the k-mer
            x[0] = (x[0] << 2 | base) & kmerMask;
            x[1] = x[1] >> 2 | (uint64_t)(3ULL - base) << bitShift;

            // once enough bases have been processed, start collecting canonical k-mers
            if (++l >= kSize)
                (x[0] <= x[1]) ? kmers.emplace_back(x[0]) : kmers.emplace_back(x[1]);
        }
        return;
    }

    // GetRCencoding will reverse complement an encoded k-mer and get it's reverse complement.
    kmer_t GetRCencoding(kmer_t encodedKmer, uint32_t kSize)
    {
        kmer_t rc = 0;
        for (uint32_t i = 0; i < kSize; i++)
        {
            rc = (rc << 2) | ((encodedKmer & 3ULL) ^ 3ULL);
            encodedKmer >>= 2;
        }
        return rc;
    }

    // DecodeKmer will decode an integer encoded k-mer.
    void DecodeKmer(kmer_t encodedKmer, uint32_t kSize, std::string& decodedKmer)
    {
        // resize the k-mer container so it can receive the decoding
        decodedKmer.resize(kSize);
        for (uint32_t i = 0; i < kSize; i++)
        {
            // add a mask to grab 2 bits at time
            uint8_t base = encodedKmer & 3ULL;
            decodedKmer[(kSize - 1) - i] = char2nt[base];
            encodedKmer >>= 2;
        }
        return;
    }

    // DecodeKmer_rc will decode an integer encoded k-mer to it's reverse complement.
    void DecodeKmer_rc(kmer_t encodedKmer, uint32_t kSize, std::string& decodedKmer)
    {
        auto rc = GetRCencoding(encodedKmer, kSize);
        DecodeKmer(rc, kSize, decodedKmer);
        return;
    }

} // namespace artic

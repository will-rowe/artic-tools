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
        uint64_t kmerMask = (0x1 << 2 * kSize) - 1;
        uint64_t bitShift = 2 * (kSize - 1);
        kmer_t fk, rk = 0;
        uint32_t l = 0;
        for (uint32_t i = 0; i < seqLen; i++)
        {

            // get 2-bit encoding of current base
            uint8_t base = nt2char[static_cast<unsigned char>(seq[i])];

            // skip non A/C/T/G bases
            if (base > 3)
                continue;

            // use the masks to add the current base to the k-mer and its reverse complement, dropping the oldest base from the k-mer
            fk = (fk << 2 | base) & kmerMask;
            rk = (rk >> 2) | (0x3 ^ base) << bitShift;

            // once enough bases have been processed, start collecting canonical k-mers
            l++;
            if (l >= kSize)
                (fk <= rk) ? kmers.push_back(fk) : kmers.push_back(rk);
        }
        return;
    }

    // GetRCencoding will reverse complement an encoded k-mer and get it's reverse complement.
    kmer_t GetRCencoding(kmer_t encodedKmer, uint32_t kSize)
    {
        kmer_t rc = 0;
        for (uint32_t i = 0; i < kSize; i++)
        {
            rc = (rc << 2) | ((encodedKmer & 0x3) ^ 0x3);
            encodedKmer >>= 2;
        }
        return rc;
    }

    // GetNextKmers will return the 4 possible next k-mers (based on actg alphabet).
    void GetNextKmers(kmer_t encodedKmer, uint32_t kSize, kmer_t& nextA, kmer_t& nextB, kmer_t& nextC, kmer_t& nextD)
    {
        uint64_t kmerMask = (0x1 << 2 * kSize) - 1;
        encodedKmer = (encodedKmer << 2) & kmerMask;
        nextA = encodedKmer;
        nextB = encodedKmer + 1;
        nextC = encodedKmer + 2;
        nextD = encodedKmer + 3;
        return;
    }

    // DecodeKmer will decode an integer encoded k-mer.
    void DecodeKmer(kmer_t encodedKmer, uint32_t kSize, std::string& decodedKmer)
    {
        // resize the k-mer container so it can receive the decoding
        decodedKmer.resize(kSize);
        for (uint32_t i = 0; i < kSize; i++)
        {
            // add a mask to grab 2 bits at time
            uint8_t base = encodedKmer & 0x3;
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

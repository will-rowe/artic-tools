#include <gtest/gtest.h>
#include <string>

#include <artic/kmers.hpp>
using namespace artic;

// some test parameters
const uint8_t kSize = 5;
const std::string sequence = "acgtana";
const std::string seqUpper = "ACGTANA";
const int expectedKmerTotal = 2; // k-mers containing N's are dropped
const kmer_t expectedEncoding = 108;
const std::string kmer1 = "ACGTA";
const std::string kmer1_rc = "TACGT";
const std::string bigSeq = "ACGTAGAAAAGG";

// reverse complementing an encoded kmer.
TEST(kmers, rc)
{
    std::string a;
    std::string b;
    artic::decodeKmer(expectedEncoding, kSize, a);
    artic::decodeKmer_rc(expectedEncoding, kSize, b);
    ASSERT_TRUE(a != b);
    ASSERT_TRUE(a == kmer1);
    ASSERT_TRUE(b == kmer1_rc);
}

// check that successive k-mer permutations returned by getNextKmers.
TEST(kmers, nextKmers)
{
    kmer_t a;
    kmer_t b;
    kmer_t c;
    kmer_t d;
    artic::getNextKmers(expectedEncoding, kSize, a, b, c, d);

    std::string as;
    std::string bs;
    std::string cs;
    std::string ds;
    artic::decodeKmer(a, kSize, as);
    artic::decodeKmer(b, kSize, bs);
    artic::decodeKmer(c, kSize, cs);
    artic::decodeKmer(d, kSize, ds);

    ASSERT_TRUE(as == ("CGTAA"));
    ASSERT_TRUE(bs == ("CGTAC"));
    ASSERT_TRUE(cs == ("CGTAG"));
    ASSERT_TRUE(ds == ("CGTAT"));
}

// encoding and decoding k-mers from a sequence.
TEST(kmers, encoding)
{
    artic::kmerset_t kmers;

    // get the first k-mer and check it produces the correct encoding
    artic::getEncodedKmers(sequence.substr(0, kSize).c_str(), kSize, kSize, kmers);
    ASSERT_TRUE(kmers.size() == 1);
    ASSERT_TRUE(kmers.find(expectedEncoding) != kmers.end());
    ASSERT_TRUE(kmers.count(expectedEncoding) == 1);

    // add the same k-mer again and check it's added to the multiset
    artic::getEncodedKmers(sequence.substr(0, kSize).c_str(), kSize, kSize, kmers);
    ASSERT_TRUE(kmers.count(expectedEncoding) == 2);

    // add the same k-mer but uppercase to check that function is case insensitive
    artic::getEncodedKmers(seqUpper.substr(0, kSize).c_str(), kSize, kSize, kmers);
    ASSERT_TRUE(kmers.count(expectedEncoding) == 3);

    // empty the set and then the whole sequence
    kmers.clear();
    artic::getEncodedKmers(sequence.c_str(), sequence.size(), kSize, kmers);
    ASSERT_TRUE(kmers.size() == expectedKmerTotal);

    // catch k-mer size error
    try
    {
        artic::getEncodedKmers(sequence.c_str(), sequence.size(), artic::MAX_K_SIZE + 1, kmers);
        FAIL() << "expected a k-mer size error";
    }
    catch (std::runtime_error& err)
    {
        std::cerr << err.what() << std::endl;
    }
    catch (...)
    {
        FAIL() << "expected a std::runtime_error";
    }

    // decode k-mers
    std::string decoded;
    std::string decoded_rc;
    for (auto kmer : kmers)
    {
        ASSERT_TRUE(kmer != 0);
        artic::decodeKmer(kmer, kSize, decoded);
        artic::decodeKmer_rc(kmer, kSize, decoded_rc);
        std::cerr << "kmer int: " << std::to_string(kmer) << " kmer strings: " << decoded << "," << decoded_rc << std::endl;
    }
}

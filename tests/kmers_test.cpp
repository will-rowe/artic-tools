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
    artic::DecodeKmer(expectedEncoding, kSize, a);
    artic::DecodeKmer_rc(expectedEncoding, kSize, b);
    ASSERT_TRUE(a != b);
    ASSERT_TRUE(a == kmer1);
    ASSERT_TRUE(b == kmer1_rc);
}

// check that successive k-mer permutations returned by GetNextKmers.
/*
TEST(kmers, nextKmers)
{
    kmer_t a;
    kmer_t b;
    kmer_t c;
    kmer_t d;
    artic::GetNextKmers(expectedEncoding, kSize, a, b, c, d);

    std::string as;
    std::string bs;
    std::string cs;
    std::string ds;
    artic::DecodeKmer(a, kSize, as);
    artic::DecodeKmer(b, kSize, bs);
    artic::DecodeKmer(c, kSize, cs);
    artic::DecodeKmer(d, kSize, ds);

    ASSERT_TRUE(as == ("CGTAA"));
    ASSERT_TRUE(bs == ("CGTAC"));
    ASSERT_TRUE(cs == ("CGTAG"));
    ASSERT_TRUE(ds == ("CGTAT"));
}
*/

// encoding and decoding k-mers from a sequence.
TEST(kmers, encoding)
{
    artic::kmerset_t kmers;

    // get the first k-mer and check it produces the correct encoding
    artic::GetEncodedKmers(sequence.substr(0, kSize).c_str(), kSize, kSize, kmers);
    ASSERT_TRUE(kmers.size() == 1);
    ASSERT_TRUE(kmers.front() == expectedEncoding);

    // check you can't reuse the provided container
    try
    {
        artic::GetEncodedKmers(sequence.substr(0, kSize).c_str(), kSize, kSize, kmers);
    }
    catch (const std::runtime_error& err)
    {
        EXPECT_EQ(err.what(), std::string("provided k-mer container already has k-mers in it"));
    }
    kmers.clear();

    // add the same k-mer but uppercase to check that function is case insensitive
    artic::GetEncodedKmers(seqUpper.substr(0, kSize).c_str(), kSize, kSize, kmers);
    ASSERT_TRUE(kmers.size() == 1);
    ASSERT_TRUE(kmers.back() == expectedEncoding);

    // empty the set and then the whole sequence
    kmers.clear();
    artic::GetEncodedKmers(sequence.c_str(), sequence.size(), kSize, kmers);
    ASSERT_TRUE(kmers.size() == expectedKmerTotal);

    // catch k-mer size error
    try
    {
        artic::GetEncodedKmers(sequence.c_str(), sequence.size(), artic::MAX_K_SIZE + 1, kmers);
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
        artic::DecodeKmer(kmer, kSize, decoded);
        artic::DecodeKmer_rc(kmer, kSize, decoded_rc);
        std::cerr << "kmer int: " << std::to_string(kmer) << " kmer strings: " << decoded << "," << decoded_rc << std::endl;
    }
}

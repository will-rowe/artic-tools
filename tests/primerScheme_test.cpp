#include <gtest/gtest.h>
#include <iostream>
#include <string>

#include <artic/log.hpp>
#include <artic/primerScheme.hpp>
using namespace artic;

// some test parameters
const unsigned int numPools = 2;
const unsigned int numPrimers = 218;
const unsigned int numAlts = 22;
const unsigned int numAmplicons = 98;
const unsigned int kSize = 17;
const std::string pool1 = "1";
const std::string pool2 = "2";
const std::string reference = std::string(TEST_DATA_PATH) + "SCoV2.reference.fasta";
const std::string refID = "MN908947.3";
const std::string inputScheme = std::string(TEST_DATA_PATH) + "SCoV2.scheme.v3.bed";

// scheme constructor
TEST(primerscheme, constructor)
{

    // catch no file
    try
    {
        auto ps = artic::PrimerScheme("");
        FAIL() << "expected a no file error";
    }
    catch (std::runtime_error& err)
    {
        EXPECT_EQ(err.what(), std::string("primer scheme input file required"));
    }
    catch (...)
    {
        FAIL() << "expected std::runtime_error";
    }

    // constructor should pass
    try
    {
        auto ps = artic::PrimerScheme(inputScheme);
    }
    catch (...)
    {
        FAIL() << "failed to construct scheme";
    }
}

// scheme validity
TEST(primerscheme, validity)
{
    artic::SchemeArgs schemeArgs;
    schemeArgs.schemeVersion = 0;
    schemeArgs.schemeFile = inputScheme;
    try
    {
        artic::Log::Init("validate_scheme");
        auto ps = artic::ValidateScheme(schemeArgs);
        EXPECT_EQ(ps.GetFileName(), inputScheme);
        EXPECT_EQ(ps.GetReferenceName(), refID);
        EXPECT_EQ(ps.GetPrimerPools().size(), numPools);
        EXPECT_EQ(ps.GetNumPrimers(), numPrimers);
        EXPECT_EQ(ps.GetNumAlts(), numAlts);
        EXPECT_EQ(ps.GetMinPrimerLen(), 22);
        EXPECT_EQ(ps.GetMaxPrimerLen(), 57);
        EXPECT_EQ(ps.GetNumAmplicons(), numAmplicons);
        EXPECT_EQ(ps.GetMeanAmpliconSpan(), 343);
    }
    catch (...)
    {
        FAIL() << "failed to construct scheme";
    }
}

// primer access
TEST(primerscheme, schemeAccess)
{
    auto ps = artic::PrimerScheme(inputScheme);

    // get a couple of primer pairs and check if they are properly paired etc.
    try
    {
        auto pp = ps.FindPrimers(40, 400);
        auto id = pp.GetName();
        ASSERT_TRUE(pp.IsProperlyPaired());
        EXPECT_EQ(id, std::string("nCoV-2019_1_LEFT_nCoV-2019_1_RIGHT"));

        // make sure pool names and IDs match up between primers and the scheme
        auto poolID = pp.GetPrimerPoolID();
        auto retPoolName = ps.GetPrimerPool(poolID);
        EXPECT_EQ(retPoolName, pool1);
        auto pp2 = ps.FindPrimers(4046, 4450);
        auto id2 = pp2.GetName();
        auto span = pp2.GetMaxSpan();
        ASSERT_TRUE(pp2.IsProperlyPaired());
        EXPECT_EQ(id2, std::string("nCoV-2019_14_LEFT_nCoV-2019_14_RIGHT"));
        EXPECT_EQ(ps.GetPrimerPool(pp2.GetPrimerPoolID()), pool2);
        EXPECT_EQ(span.first, 4044);
        EXPECT_EQ(span.second, 4450);
    }
    catch (std::runtime_error& err)
    {
        FAIL() << "runtime error: " << err.what();
    }
    try
    {
        auto pp = ps.FindPrimers(300, 400);
        ASSERT_FALSE(pp.IsProperlyPaired());
        EXPECT_EQ(ps.GetPrimerPool(pp.GetPrimerPoolID()), std::string("unmatched"));
    }
    catch (std::runtime_error& err)
    {
        FAIL() << "runtime error: " << err.what();
    }
    auto primerPools = ps.GetPrimerPools();
    EXPECT_EQ(primerPools.size(), numPools);
}

// primer sites
TEST(primerscheme, primerSites)
{
    auto ps = artic::PrimerScheme(inputScheme);
    bool primerSite;

    // check out of bounds
    try
    {
        primerSite = ps.CheckPrimerSite(0);
    }
    catch (std::runtime_error& err)
    {
        EXPECT_STREQ("query position outside of primer scheme bounds", err.what());
    }

    // check all primers are returned as being in the scheme, and check that an offset falls outside
    auto amplicons = ps.GetExpAmplicons();
    for (auto amplicon : amplicons)
    {

        // check the forward and reverse primers are in primer sites register
        ASSERT_TRUE(ps.CheckPrimerSite(amplicon.GetForwardPrimer()->GetStart()));
        ASSERT_TRUE(ps.CheckPrimerSite(amplicon.GetReversePrimer()->GetEnd() - 1));

        // check the just outside the primer is not in the register
        ASSERT_FALSE(ps.CheckPrimerSite(amplicon.GetForwardPrimer()->GetEnd()));
        ASSERT_FALSE(ps.CheckPrimerSite(amplicon.GetReversePrimer()->GetStart() - 1));
    }
}

// primer sequence
TEST(primerscheme, primerSeq)
{
    auto ps = artic::PrimerScheme(inputScheme);
    auto pp = ps.FindPrimers(40, 400);
    auto p1 = pp.GetForwardPrimer();
    faidx_t* fai = fai_load(reference.c_str());
    std::string seq;
    p1->GetSeq(fai, ps.GetReferenceName(), seq);
    if (fai)
        fai_destroy(fai);
    EXPECT_EQ(seq.size(), p1->GetLen());
    EXPECT_STREQ("ACCAACCAACTTTCGATCTCTTGT", seq.c_str());
}

// expected amplicons
TEST(primerscheme, amplicons)
{
    auto ps = artic::PrimerScheme(inputScheme);
    auto amplicons = ps.GetExpAmplicons();
    EXPECT_EQ(amplicons.size(), ps.GetNumAmplicons());
}

// primer kmers
TEST(primerscheme, kmers)
{

    // load a scheme, check for k-mer size, populate k-mer map
    auto ps = artic::PrimerScheme(inputScheme);
    artic::kmermap_t kmerMap;
    try
    {
        ps.GetPrimerKmers(reference, 100, kmerMap);
        FAIL() << "expected a k-mer size error";
    }
    catch (std::runtime_error& err)
    {
        EXPECT_EQ(err.what(), std::string("requested k-mer size is greater than the shortest primer in the scheme (22)"));
    }

    ps.GetPrimerKmers(reference, kSize, kmerMap);

    // load the another copy of the scheme and check primer k-mer lookup
    auto ps2 = artic::PrimerScheme(inputScheme);
    EXPECT_EQ(ps.GetNumAmplicons(), ps2.GetNumAmplicons());
    faidx_t* fai = fai_load(reference.c_str());
    std::string seq1;
    std::string seq2;
    artic::kmerset_t kmers;
    for (auto amplicon : ps2.GetExpAmplicons())
    {

        // get k-mers from forward and reverse primers
        auto p1 = amplicon.GetForwardPrimer();
        p1->GetSeq(fai, ps2.GetReferenceName(), seq1);
        artic::GetEncodedKmers(seq1.c_str(), seq1.size(), kSize, kmers);
        EXPECT_EQ(kmers.size(), (seq1.size() - kSize + 1));
        auto p2 = amplicon.GetReversePrimer();
        p2->GetSeq(fai, ps2.GetReferenceName(), seq2);
        artic::GetEncodedKmers(seq2.c_str(), seq2.size(), kSize, kmers);
        EXPECT_EQ(kmers.size(), ((seq1.size() - kSize + 1) + (seq2.size() - kSize + 1)));

        // loop over the amplicon primer k-mers and check against the other scheme copy
        for (auto kmer : kmers)
        {
            auto it = kmerMap.find(kmer);
            ASSERT_TRUE(it != kmerMap.end());
            bool found = false;
            for (auto id : it->second)
            {
                if (id == amplicon.GetID())
                {
                    found = true;
                    EXPECT_EQ(ps.GetAmpliconName(id), amplicon.GetName());
                }
            }
            ASSERT_TRUE(found);
        }
        kmers.clear();
        seq1.clear();
        seq2.clear();
    }
    if (fai)
        fai_destroy(fai);

    /*
    for (auto x : kmerMap)
    {
        std::string decoded;
        artic::DecodeKmer(x.first, kSize, decoded);
        std::cerr << "kmer hash: " << x.first << " --- " << decoded << std::endl;
        for (auto y : x.second)
        {
            std::cerr << "\tlocated in amplicon: " << y << " --- " << ps.GetAmpliconName(y) << std::endl;
        }
    }
    */
}

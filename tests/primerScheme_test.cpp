#include <gtest/gtest.h>
#include <iostream>
#include <string>

#include <artic/primerScheme.hpp>
using namespace artic;

// some test parameters
const unsigned int version = 3;
const unsigned int numPools = 2;
const unsigned int numPrimers = 218;
const unsigned int numAlts = 22;
const unsigned int numAmplicons = 98;
const std::string pool1 = "nCoV-2019_1";
const std::string pool2 = "nCoV-2019_2";

const std::string inputScheme = std::string(TEST_DATA_PATH) + "SCoV2.scheme.v3.bed";

// scheme constructor
TEST(primerscheme, constructor)
{

    // catch no file
    try
    {
        artic::PrimerScheme("", version);
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

    // catch bad version
    try
    {
        artic::PrimerScheme(inputScheme, 666);
        FAIL() << "expected a no file error";
    }
    catch (std::runtime_error& err)
    {
        EXPECT_EQ(err.what(), std::string("unrecognised primer scheme version: 666"));
    }
    catch (...)
    {
        FAIL() << "expected std::runtime_error";
    }

    // constructor should pass
    try
    {
        artic::PrimerScheme(inputScheme, 3);
    }
    catch (...)
    {
        FAIL() << "failed to construct scheme";
    }
}

// scheme validity
TEST(primerscheme, validity)
{
    auto ps = artic::PrimerScheme(inputScheme, 3);
    EXPECT_EQ(ps.GetVersion(), version);
    EXPECT_EQ(ps.GetNumPrimers(), numPrimers);
    EXPECT_EQ(ps.GetNumAlts(), numAlts);
    EXPECT_EQ(ps.GetNumAmplicons(), numAmplicons);
    EXPECT_EQ(ps.GetPrimerPools().size(), numPools);
}

// primer access
TEST(primerscheme, access)
{
    auto ps = artic::PrimerScheme(inputScheme, 3);

    // get a couple of primer pairs and check if they are properly paired etc.
    try
    {
        auto pp = ps.FindPrimers(40, 400);
        auto id = pp.GetID();
        ASSERT_TRUE(pp.IsProperlyPaired());
        EXPECT_EQ(id, std::string("nCoV-2019_1_LEFT_nCoV-2019_1_RIGHT"));

        // make sure pool names and IDs match up between primers and the scheme
        auto poolName = pp.GetPrimerPool();
        auto retPoolID = ps.GetPrimerPoolID(poolName);
        auto retPoolName = ps.GetPrimerPool(retPoolID);
        EXPECT_EQ(pp.GetPrimerPool(), pool1);
        EXPECT_EQ(retPoolName, pool1);

        auto pp2 = ps.FindPrimers(4046, 4450);
        auto id2 = pp2.GetID();
        auto span = pp2.GetMaxSpan();
        std::cout << id2 << std::endl;
        ASSERT_TRUE(pp2.IsProperlyPaired());
        EXPECT_EQ(id2, std::string("nCoV-2019_14_LEFT_nCoV-2019_14_RIGHT"));
        EXPECT_EQ(pp2.GetPrimerPool(), pool2);
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
        EXPECT_EQ(pp.GetPrimerPool(), std::string("unmatched"));
    }
    catch (std::runtime_error& err)
    {
        FAIL() << "runtime error: " << err.what();
    }
    auto primerPools = ps.GetPrimerPools();
    EXPECT_EQ(primerPools.size(), 2);
}

// primer sites
TEST(primerscheme, primersites)
{
    auto ps = artic::PrimerScheme(inputScheme, 3);
    bool primerSite;

    // check out of bounds
    try
    {
        primerSite = ps.CheckPrimerSite(0, pool1);
    }
    catch (std::runtime_error& err)
    {
        EXPECT_STREQ("query position outside of primer scheme bounds", err.what());
    }

    // basic check for primer sites
    primerSite = ps.CheckPrimerSite(40, pool1);
    ASSERT_TRUE(primerSite);
    primerSite = ps.CheckPrimerSite(40, pool2);
    ASSERT_FALSE(primerSite);
}
#include <gtest/gtest.h>
#include <iostream>
#include <string>

#include <artic/log.hpp>
#include <artic/vcfCheck.hpp>
using namespace artic;

// some test parameters
const std::string inputScheme = std::string(TEST_DATA_PATH) + "SCoV2.scheme.v3.bed";
const std::string vcfIn = std::string(TEST_DATA_PATH) + "CVR1.merged.vcf.gz";
const std::string reportOut = std::string(TEST_DATA_PATH) + "CVR1.merged.vcf.report.txt";
const unsigned int minVarQual = 1;

// some expected counts
const unsigned int totalVars = 22;
const unsigned int primerSeqVars = 1;
const unsigned int overlapVars = 5;
const unsigned int lowQualVars = 0;

// vcfChecker constructor
TEST(vcfchecker, run)
{
    try
    {
        artic::Log::Init("check_vcf");
        auto ps = artic::PrimerScheme(inputScheme);
        auto checker = artic::VcfChecker(&ps, vcfIn, reportOut, "", minVarQual);
        checker.Run();
        EXPECT_EQ(totalVars, checker.GetNumRecords());
        EXPECT_EQ(primerSeqVars, checker.GetNumInPrimerSite());
        EXPECT_EQ(overlapVars, checker.GetNumInOverlap());
        EXPECT_EQ(lowQualVars, checker.GetNumLowQual());
    }
    catch (...)
    {
        FAIL() << "failed to check vcf";
    }
}

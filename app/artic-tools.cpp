#include <CLI/CLI.hpp>
#include <htslib/faidx.h>
#include <string>
#include <vector>

#include <artic/amplitig.hpp>
#include <artic/log.hpp>
#include <artic/primerScheme.hpp>
#include <artic/softmask.hpp>
#include <artic/vcfCheck.hpp>
#include <artic/version.hpp>
using namespace artic;

int main(int argc, char** argv)
{

    // get the user command as a string for logging and BAM headers
    std::stringstream userCmd;
    for (int i = 1; i < argc; ++i)
    {
        if (i != 1)
            userCmd << " ";
        userCmd << argv[i];
    }

    // set up the app CLI
    CLI::App app{PROG_NAME + " is a set of artic pipeline utilities"};
    auto versionCallback = [](int) { std::cout << artic::GetVersion() << std::endl; exit(0); };
    app.add_flag_function("-v,--version", versionCallback, "Print version and exit");

    // set up the subcommands
    app.require_subcommand();
    CLI::App* softmaskCmd = app.add_subcommand("align_trim", "Trim alignments from an amplicon scheme");
    CLI::App* getterCmd = app.add_subcommand("get_scheme", "Download an ARTIC primer scheme and reference sequence");
    CLI::App* validatorCmd = app.add_subcommand("validate_scheme", "Validate an amplicon scheme for compliance with ARTIC standards");
    CLI::App* vcfFilterCmd = app.add_subcommand("check_vcf", "Check a VCF file based on primer scheme info and user-defined cut offs");
    CLI::App* amplitigCmd = app.add_subcommand("get_amplitigs", "Generate amplitigs from a reference alignment");

    // set up a struct to pass arguments
    // TODO: have a constructor for defaults?
    artic::SchemeArgs schemeArgs;
    schemeArgs.schemeVersion = 0; // indicates latest scheme available

    // add softmask/amplitiger options and flags
    std::string inputFile;
    std::vector<std::string> inputFiles;
    std::string outFileName;
    unsigned int minMAPQ = 15;
    unsigned int normalise = 100;
    unsigned int kmerSize = 21;
    bool primerStart = false;
    bool removeBadPairs = false;
    bool noReadGroups = false;
    bool verbose = false;
    softmaskCmd->add_option("-b,--inputFile", inputFile, "The input BAM file (will try STDIN if not provided)");
    softmaskCmd->add_option("scheme", schemeArgs.schemeFile, "The ARTIC primer scheme")->required()->check(CLI::ExistingFile);
    softmaskCmd->add_option("--minMAPQ", minMAPQ, "A minimum MAPQ threshold for processing alignments (default = 15)");
    softmaskCmd->add_option("--normalise", normalise, "Subsample to N coverage per strand (default = 100, deactivate with 0)");
    softmaskCmd->add_option("--report", outFileName, "Output an align_trim report to file");
    softmaskCmd->add_flag("--start", primerStart, "Trim to start of primers instead of ends");
    softmaskCmd->add_flag("--remove-incorrect-pairs", removeBadPairs, "Remove amplicons with incorrect primer pairs");
    softmaskCmd->add_flag("--no-read-groups", noReadGroups, "Do not divide reads into groups in SAM output");
    softmaskCmd->add_flag("--verbose", verbose, "Output debugging information to STDERR");

    // add amplitig options and flags
    amplitigCmd->add_option("-i,--fastqFiles", inputFiles, "The input FASTQ files");
    amplitigCmd->add_option("scheme", schemeArgs.schemeFile, "The ARTIC primer scheme")->required()->check(CLI::ExistingFile);
    amplitigCmd->add_option("-r,--refSeq", schemeArgs.refSeqFile, "The reference sequence for the primer scheme (FASTA format)")->required();
    amplitigCmd->add_option("-k,--kmerSize", kmerSize, "The k-mer size to use (default = 21)");
    amplitigCmd->add_flag("--verbose", verbose, "Output debugging information to STDERR");

    // add get options and flags
    getterCmd->add_option("scheme", schemeArgs.schemeName, "The name of the scheme to download (ebola|nipah|scov2)")->required();
    getterCmd->add_option("--schemeVersion", schemeArgs.schemeVersion, "The ARTIC primer scheme version (default = latest)")->default_val(0);
    getterCmd->add_option("-o,--outDir", schemeArgs.outDir, "The directory to write the scheme and reference sequence to");

    // add validator options and flags
    validatorCmd->add_option("scheme", schemeArgs.schemeFile, "The primer scheme to validate")->required()->check(CLI::ExistingFile);
    validatorCmd->add_option("--schemeVersion", schemeArgs.schemeVersion, "The ARTIC primer scheme version (default = latest)");
    validatorCmd->add_option("-o,--outputPrimerSeqs", schemeArgs.primerSeqsFile, "If provided, will write primer sequences as multiFASTA (requires --refSeq to be provided)");
    validatorCmd->add_option("-r,--refSeq", schemeArgs.refSeqFile, "The reference sequence for the primer scheme (FASTA format)");
    validatorCmd->add_option("--outputInserts", schemeArgs.insertsFile, "If provided, will write primer scheme inserts as BED (exluding primer sequences)");

    // add vcfFilter options and flags
    std::string vcfIn;
    bool dropPrimerVars = false;
    bool dropOverlapFails = false;
    vcfFilterCmd->add_option("vcf", vcfIn, "The input VCF file to filter")->required()->check(CLI::ExistingFile);
    vcfFilterCmd->add_option("scheme", schemeArgs.schemeName, "The primer scheme to use")->required()->check(CLI::ExistingFile);
    vcfFilterCmd->add_option("-o,--vcfOut", outFileName, "If provided, will write variants that pass checks");
    vcfFilterCmd->add_flag("--dropPrimerVars", dropPrimerVars, "Will drop variants called within primer regions for the pool");
    vcfFilterCmd->add_flag("--dropOverlapFails", dropOverlapFails, "Will drop variants called once within amplicon overlap regions");

    // add the align_trim callback
    // 1. validate the primer scheme
    // 2. run the softmasker
    softmaskCmd->callback([&]() {
        artic::Log::Init("align_trim");
        LOG_TRACE("starting align trim");
        auto ps = artic::ValidateScheme(schemeArgs);
        auto masker = artic::Softmasker(&ps, inputFile, userCmd.str(), minMAPQ, normalise, removeBadPairs, noReadGroups, primerStart, outFileName);
        masker.Run(verbose);
    });

    // add the amplitigger callback
    // 1. vadlidate the primer scheme
    // 2. collect the primer k-mers
    // 2. run the amplitiger
    amplitigCmd->callback([&]() {
        artic::Log::Init("get_amplitigs");
        LOG_TRACE("starting amplitigger");
        auto ps = artic::ValidateScheme(schemeArgs);
        auto amplitigger = artic::Amplitigger(&ps, schemeArgs.refSeqFile, inputFiles, userCmd.str(), kmerSize);
        amplitigger.Run(verbose);
    });

    // add the getter callback
    // 1. download the scheme
    // 2. validate the scheme
    getterCmd->callback([&]() {
        artic::Log::Init("get_scheme");
        LOG_TRACE("starting primer scheme downloader");
        artic::DownloadScheme(schemeArgs);
        artic::ValidateScheme(schemeArgs);
    });

    // add the validator callback
    // 1. validate the scheme
    validatorCmd->callback([&]() {
        artic::Log::Init("validate_scheme");
        LOG_TRACE("starting primer scheme validator");
        artic::ValidateScheme(schemeArgs);
    });

    // add the vcfFilter callback
    // 1. run the vcf filterer
    vcfFilterCmd->callback([&]() {
        artic::Log::Init("check_vcf");
        LOG_TRACE("starting VCF checker");
        auto ps = artic::ValidateScheme(schemeArgs);
        auto filter = artic::VcfChecker(&ps, vcfIn, outFileName, dropPrimerVars, dropOverlapFails);
        filter.Run();
    });

    // parse CLI
    try
    {
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError& e)
    {
        return app.exit(e);
    }
    catch (const std::runtime_error& re)
    {
        std::cerr << "error--> " << re.what() << std::endl;
        return -1;
    }
    catch (...)
    {
        std::cerr << "unknown error" << std::endl;
        return -1;
    }
    return 0;
}

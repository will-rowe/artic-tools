#include <CLI/CLI.hpp>
#include <string>
#include <vector>

#include <artic/primerScheme.hpp>
#include <artic/softmask.hpp>
#include <artic/version.hpp>
using namespace artic;

int main(int argc, char** argv)
{

    // set up the app CLI
    CLI::App app{PROG_NAME + " is a set of artic pipeline utilities"};
    auto versionCallback = [](int) { std::cout << artic::GetVersion() << std::endl; exit(0); };
    app.add_flag_function("-v,--version", versionCallback, "Print version and exit");

    // set up the subcommands
    app.require_subcommand();
    int errorCode = 0;

    // add softmask subcmd
    CLI::App* softmaskCmd = app.add_subcommand("align_trim", "Trim alignments from an amplicon scheme");

    // add softmask required fields
    std::string bamFile;
    std::string primerSchemeFile;
    softmaskCmd->add_option("-b,--bamFile", bamFile, "The input bam file (will try STDIN if not provided)");
    softmaskCmd->add_option("-s,--primerScheme", primerSchemeFile, "The ARTIC primer scheme")->required();

    // add softmask optional fields
    int primerSchemeVersion = 3;
    unsigned int minMAPQ = 15;
    unsigned int normalise = 100;
    std::string report;
    bool primerStart = false;
    bool removeBadPairs = false;
    bool noReadGroups = false;
    bool verbose = false;
    softmaskCmd->add_option("--primerSchemeVersion", primerSchemeVersion, "The ARTIC primer scheme version (default = 3)");
    softmaskCmd->add_option("--minMAPQ", minMAPQ, "A minimum MAPQ threshold for processing alignments (default = 15)");
    softmaskCmd->add_option("--normalise", normalise, "Subsample to N coverage per strand (default = 100, deactivate with 0)");
    softmaskCmd->add_option("--report", report, "Output an align_trim report to file");
    softmaskCmd->add_flag("--start", primerStart, "Trim to start of primers instead of ends");
    softmaskCmd->add_flag("--remove-incorrect-pairs", removeBadPairs, "Remove amplicons with incorrect primer pairs");
    softmaskCmd->add_flag("--no-read-groups", noReadGroups, "Do not divide reads into groups in SAM output");
    softmaskCmd->add_flag("--verbose", verbose, "Output debugging information to STDERR");

    // get the user command as a string for logging and BAM headers
    std::stringstream userCmd;
    for (int i = 1; i < argc; ++i)
    {
        if (i != 1)
            userCmd << " ";
        userCmd << argv[i];
    }

    // add the softmask callback
    softmaskCmd->callback([&]() {
        try
        {
            // load and check the primer scheme
            artic::PrimerScheme primerScheme = artic::PrimerScheme(primerSchemeFile, primerSchemeVersion);

            // setup and run the softmasker
            auto masker = artic::Softmasker(&primerScheme, bamFile, userCmd.str(), minMAPQ, normalise, removeBadPairs, noReadGroups, primerStart, report);
            masker.Run(verbose);
        }
        catch (const std::runtime_error& re)
        {
            std::cerr << re.what() << std::endl;
            errorCode = -1;
        }
        catch (const std::exception& ex)
        {
            std::cerr << "exception: " << ex.what() << std::endl;
            errorCode = -1;
        }
        catch (...)
        {
            std::cerr << "unknown error" << std::endl;
            errorCode = -1;
        }
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

    return errorCode;
}

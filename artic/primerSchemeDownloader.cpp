#include <boost/filesystem.hpp>
#include <iostream>
#include <map>
#include <sstream>

#include "log.hpp"
#include "primerScheme.hpp"

// Max version available for each ARTIC primer scheme
const uint maxEbolaVersion = 3;
const uint maxNipahVersion = 1;
const uint maxSARSCoV2Version = 3;

// File helpers
const std::string schemeExt = ".primer.bed";
const std::string refExt = ".reference.fasta";
const std::string baseURL = "https://github.com/artic-network/primer-schemes/raw/master/";
const std::pair<std::string, std::string> ebovTags = {"ZaireEbola", "ZaireEbola"};
const std::pair<std::string, std::string> nivTags = {"Nipah", "NiV_6_Malaysia"};
const std::pair<std::string, std::string> scovTags = {"nCoV-2019", "nCoV-2019"};

// Scheme enum is used in the switch statement which constructs the download link
enum Scheme
{
    InvalidScheme,
    Ebov,
    Nipah,
    Scov2,
};

Scheme resolveScheme(std::string query)
{
    static const std::map<std::string, Scheme> schemeStrings{
        {"ebola", Ebov},
        {"nipah", Nipah},
        {"scov2", Scov2},
        {"ncov", Scov2}, // added for backward compatability
    };
    auto itr = schemeStrings.find(query);
    if (itr != schemeStrings.end())
        return itr->second;
    return InvalidScheme;
}

std::string resolveLink(const std::string& tag1, unsigned int version, const std::string& tag2)
{
    std::stringstream s;
    s << baseURL << tag1 << "/V" << version << "/" << tag2;
    return s.str();
}

// DownloadScheme will download a specified primer scheme and the reference sequence.
// if 0 is provided as a version, or an unknown version provided, the latest scheme version will be used.
void artic::DownloadScheme(SchemeArgs& args)
{
    std::string schemeURL;
    unsigned int version = args.schemeVersion;
    switch (resolveScheme(args.schemeName))
    {
        case Ebov:
            if (version == 0 || version > maxEbolaVersion)
                version = maxEbolaVersion;
            schemeURL = resolveLink(ebovTags.first, version, ebovTags.second);
            break;
        case Nipah:
            if (version == 0 || version > maxNipahVersion)
                version = maxNipahVersion;
            schemeURL = resolveLink(nivTags.first, version, nivTags.second);
            break;
        case Scov2:
            if (version == 0 || version > maxSARSCoV2Version)
                version = maxSARSCoV2Version;
            schemeURL = resolveLink(scovTags.first, version, scovTags.second);
            break;
        default:
            throw std::runtime_error("unknown scheme: " + args.schemeName);
    }
    LOG_TRACE("\tscheme: {}", args.schemeName);
    LOG_TRACE("\tversion: {}", args.schemeVersion);
    if ((args.schemeVersion != version) && (args.schemeVersion != 0))
        LOG_WARN("requested scheme version not found, using latest version instead (v{})", version);
    if (args.schemeVersion == 0)
        LOG_WARN("requested latest scheme version, using v{}", version);
    args.schemeVersion = version;

    // create filenames
    std::stringstream of;
    if (args.outDir.size() != 0)
    {
        if (!boost::filesystem::is_directory(args.outDir) || !boost::filesystem::exists(args.outDir))
            boost::filesystem::create_directories(args.outDir);
        of << args.outDir << "/";
    }
    of << args.schemeName << ".v" << args.schemeVersion;

    // download using wget sys call
    LOG_TRACE("downloading data")
    std::stringstream baseCmd;
    baseCmd << "wget -q -O " << args.schemeName << ".v" << args.schemeVersion;
    std::string cmd1 = "wget -q -O " + of.str() + schemeExt + " " + schemeURL + schemeExt + " 2>/dev/null";
    std::string cmd2 = "wget -q -O " + of.str() + refExt + " " + schemeURL + refExt + " 2>/dev/null";
    int res;
    res = system(cmd1.c_str());
    if (res != EXIT_SUCCESS)
        throw std::runtime_error("could not download scheme with wget");
    res = system(cmd2.c_str());
    if (res != EXIT_SUCCESS)
        throw std::runtime_error("could not download reference with wget");
    LOG_TRACE("\t{}{}", of.str(), schemeExt);
    LOG_TRACE("\t{}{}", of.str(), refExt);

    // pass back the downloaded scheme via the args struct
    args.schemeFile = of.str() + schemeExt;
    args.refSeqFile = of.str() + refExt;
}

#include <boost/filesystem.hpp>
#include <iostream>
#include <map>

#include "log.hpp"
#include "primerScheme.hpp"

// ARTIC scheme downloads
// Note: it would be nice to have all these in one place and under a unified naming system...
const std::vector<std::string> ebovUrls = {"https://github.com/artic-network/primer-schemes/raw/master/ZaireEbola/V1/ZaireEbola.scheme.bed", "https://github.com/artic-network/primer-schemes/raw/master/ZaireEbola/V2/ZaireEbola.scheme.bed", "https://github.com/artic-network/primer-schemes/raw/master/ZaireEbola/V3/Ebov-10-Pan.bed"};
const std::vector<std::string> nivUrls = {"https://github.com/artic-network/primer-schemes/raw/master/Nipah/V1/NiV_6_Malaysia.bed"};
const std::vector<std::string> scovUrls = {"https://github.com/artic-network/artic-ncov2019/raw/master/primer_schemes/nCoV-2019/V1/nCoV-2019.scheme.bed", "https://github.com/artic-network/artic-ncov2019/raw/master/primer_schemes/nCoV-2019/V2/nCoV-2019.scheme.bed", "https://github.com/artic-network/artic-ncov2019/raw/master/primer_schemes/nCoV-2019/V3/nCoV-2019.scheme.bed"};

const std::vector<std::string> ebovRefUrls = {"https://github.com/artic-network/primer-schemes/raw/master/ZaireEbola/V1/ZaireEbola.reference.fasta", "https://github.com/artic-network/primer-schemes/raw/master/ZaireEbola/V2/ZaireEbola.reference.fasta", "https://github.com/artic-network/primer-schemes/raw/master/ZaireEbola/V3/Ebov-10-Pan.fasta"};
const std::vector<std::string> nivRefUrls = {"https://github.com/artic-network/primer-schemes/raw/master/Nipah/V1/NiV_6_Malaysia.fasta"};
const std::vector<std::string> scovRefUrls = {"https://github.com/artic-network/artic-ncov2019/raw/master/primer_schemes/nCoV-2019/V1/nCoV-2019.reference.fasta", "https://github.com/artic-network/artic-ncov2019/raw/master/primer_schemes/nCoV-2019/V2/nCoV-2019.reference.fasta", "https://github.com/artic-network/artic-ncov2019/raw/master/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta"};

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
    };
    auto itr = schemeStrings.find(query);
    if (itr != schemeStrings.end())
        return itr->second;
    return InvalidScheme;
}

// DownloadScheme will download a specified primer scheme and the reference sequence.
// if 0 is provided as a version, or an unknown version provided, the latest scheme version will be used.
void artic::DownloadScheme(const std::string& schemeName, unsigned int requestedVersion, const std::string& outDir)
{
    std::string schemeURL;
    std::string refURL;
    unsigned int version = requestedVersion;
    switch (resolveScheme(schemeName))
    {
        case Ebov:
            if (version == 0 || version > ebovUrls.size())
                version = ebovUrls.size();
            schemeURL = ebovUrls.at(version - 1);
            refURL = ebovRefUrls.at(version - 1);
            break;
        case Nipah:
            if (version == 0 || version > nivUrls.size())
                version = nivUrls.size();
            schemeURL = nivUrls.at(version - 1);
            refURL = nivRefUrls.at(version - 1);
            break;
        case Scov2:
            if (version == 0 || version > scovUrls.size())
                version = scovUrls.size();
            schemeURL = scovUrls.at(version - 1);
            refURL = scovRefUrls.at(version - 1);
            break;
        default:
            throw std::runtime_error("unknown scheme: " + schemeName);
    }
    artic::Log::Init("getter");
    LOG_TRACE("starting scheme getter");
    LOG_TRACE("\tscheme: {}", schemeName);
    LOG_TRACE("\tversion: {}", requestedVersion);
    if ((requestedVersion != version) && (requestedVersion != 0))
        LOG_WARN("requested scheme version not found, using latest version instead (v{})", version);

    // create filenames
    std::stringstream of1;
    std::stringstream of2;
    if (outDir.size() != 0)
    {
        if (!boost::filesystem::is_directory(outDir) || !boost::filesystem::exists(outDir))
            boost::filesystem::create_directory(outDir);
        of1 << outDir << "/";
        of2 << outDir << "/";
    }
    of1 << schemeName << ".v" << version << ".scheme.bed";
    of2 << schemeName << ".v" << version << ".reference.fasta";

    // download using wget sys call
    LOG_TRACE("downloading data")
    std::stringstream baseCmd;
    baseCmd << "wget -q -O ";
    baseCmd << schemeName << ".v" << version;
    std::string cmd1 = "wget -q -O " + of1.str() + " " + schemeURL + " 2>/dev/null";
    std::string cmd2 = "wget -q -O " + of2.str() + " " + refURL + " 2>/dev/null";
    int res;
    res = system(cmd1.c_str());
    if (res != EXIT_SUCCESS)
        throw std::runtime_error("could not download scheme with wget");
    res = system(cmd2.c_str());
    if (res != EXIT_SUCCESS)
        throw std::runtime_error("could not download reference with wget");
    LOG_TRACE("\t{}", of1.str());
    LOG_TRACE("\t{}", of2.str());

    // validate the scheme
    LOG_TRACE("checking scheme")
    auto ps = artic::PrimerScheme(of1.str(), version);
}

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <curl/curl.h>
#include <iostream>
#include <picosha2.h>
#include <sstream>

#include "log.hpp"
#include "primerScheme.hpp"

namespace pt = boost::property_tree;

// ARTIC_MANIFEST_URL contains all the metadata for the ARTIC Network primer schemes.
const std::string ARTIC_MANIFEST_URL = "https://raw.githubusercontent.com/artic-network/primer-schemes/master/schemes_manifest.json";
const std::string SCHEME_EXT = ".primer.bed";
const std::string REF_EXT = ".reference.fasta";

// write2stream is used by curl to download to a stream.
size_t write2stream(char* ptr, size_t size, size_t nmemb, void* userdata)
{
    std::stringstream* stream = (std::stringstream*)userdata;
    size_t count = size * nmemb;
    stream->write(ptr, count);
    return count;
}

// write2disk is used by curl to download to disk.
size_t write2disk(void* ptr, size_t size, size_t nmemb, FILE* stream)
{
    size_t count = fwrite(ptr, size, nmemb, stream);
    return count;
}

// downloadJSON will download a JSON file from a URL and load it as a JSON object.
void downloadJSON(const std::string& url, pt::ptree& jsonTree)
{

    // set up curl to write to a stream
    CURL* curl = curl_easy_init();
    if (!curl)
        throw std::runtime_error("could not init curl for file download");
    std::stringstream stream;
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write2stream);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &stream);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_MAXREDIRS, 3L);

    // download the manifest to stream
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    auto res = curl_easy_perform(curl);
    if (res != 0)
        throw std::runtime_error("could not download ARTIC manifest: " + std::to_string(res));
    curl_easy_cleanup(curl);

    // deserialize from curl stream
    pt::read_json(stream, jsonTree);
    return;
}

// downloadFile will download a file from a URL and save it on disk.
void downloadFile(const std::string& url, const std::string& fname)
{
    CURL* curl = curl_easy_init();
    if (!curl)
        throw std::runtime_error("could not init curl for file download");
    FILE* fp = fopen(fname.c_str(), "wb");
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write2disk);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_MAXREDIRS, 3L);
    auto res = curl_easy_perform(curl);
    if (res != 0)
        throw std::runtime_error("could not download file: " + url + " (" + std::to_string(res) + ")");
    curl_easy_cleanup(curl);
    fclose(fp);
}

// DownloadScheme will download a specified primer scheme and the reference sequence.
// if 0 is provided as a version, or an unknown version provided, the latest scheme version will be used.
void artic::DownloadScheme(SchemeArgs& args)
{
    LOG_TRACE("\trequested scheme:\t{}", args.schemeName);
    if (args.schemeVersion == 0)
    {
        LOG_TRACE("\trequested version:\tlatest");
    }
    else
    {
        LOG_TRACE("\trequested version:\t{}", args.schemeVersion);
    }
    LOG_TRACE("fetching manifest file");
    LOG_TRACE("\tARTIC manifest URL:\t{}", ARTIC_MANIFEST_URL);
    pt::ptree manifest;
    downloadJSON(ARTIC_MANIFEST_URL, manifest);
    LOG_TRACE("\tARTIC repository DOI:\t{}", manifest.get<std::string>("latest_doi"));

    // grab the schemes as a subtree and then check each for the requested scheme (via the alias field)
    LOG_TRACE("finding primer scheme");
    std::vector<std::string> checkedAliases;
    std::pair<std::string, std::string> primer_url;
    std::pair<std::string, std::string> reference_url;
    pt::ptree& schemes = manifest.get_child("schemes");
    for (const auto& scheme : schemes)
    {
        for (const auto& alias : scheme.second.get_child("aliases"))
        {
            checkedAliases.push_back(alias.second.data());
            if (!boost::iequals(checkedAliases.back(), args.schemeName))
                continue;

            // check available versions and get the primer_url
            LOG_TRACE("\tfound requested scheme:\t{} (using alias {})", scheme.first, args.schemeName);
            unsigned int latestVersion = scheme.second.get<unsigned int>("latest_version");
            if (args.schemeVersion > latestVersion)
            {
                LOG_WARN("\trequested version not found (v{}), using latest version instead (v{})", args.schemeVersion, latestVersion);
                args.schemeVersion = latestVersion;
            }
            if (args.schemeVersion == 0)
                args.schemeVersion = latestVersion;
            primer_url.first = scheme.second.get_child("primer_urls").get<std::string>(std::to_string(args.schemeVersion));
            primer_url.second = scheme.second.get_child("primer_sha256_checksums").get<std::string>(std::to_string(args.schemeVersion));
            reference_url.first = scheme.second.get_child("reference_urls").get<std::string>(std::to_string(args.schemeVersion));
            reference_url.second = scheme.second.get_child("reference_sha256_checksums").get<std::string>(std::to_string(args.schemeVersion));
        }
    }

    // check that a scheme was found
    if (primer_url.first.empty() || reference_url.first.empty())
    {
        LOG_WARN("\tscheme not found:\t{}", args.schemeName);
        LOG_WARN("listing available scheme aliases (case insensitive)", args.schemeName);
        for (const auto& i : checkedAliases)
            LOG_WARN("\t- {}", i);
        throw std::runtime_error("no primer scheme available for " + args.schemeName);
    }

    // create filenames
    std::stringstream fp;
    if (args.outDir.size() != 0)
    {
        if (!boost::filesystem::is_directory(args.outDir) || !boost::filesystem::exists(args.outDir))
            boost::filesystem::create_directories(args.outDir);
        fp << args.outDir << "/";
    }
    fp << args.schemeName << ".v" << args.schemeVersion;
    auto f1 = fp.str() + SCHEME_EXT;
    auto f2 = fp.str() + REF_EXT;

    // download files
    LOG_TRACE("downloading primer scheme")
    downloadFile(primer_url.first, f1);
    downloadFile(reference_url.first, f2);
    LOG_TRACE("\tsaving primers to:\t{}", f1);
    LOG_TRACE("\tsaving reference to:\t{}", f2);

    // check integrity
    LOG_TRACE("comparing checksums")
    std::ifstream f1Check(f1, std::ios::binary);
    std::vector<unsigned char> s1(picosha2::k_digest_size);
    picosha2::hash256(f1Check, s1.begin(), s1.end());
    std::string f1h;
    picosha2::bytes_to_hex_string(s1.begin(), s1.end(), f1h);
    LOG_TRACE("\tsha256 for primers:\t{}", f1h);
    if (f1h != primer_url.second)
        throw std::runtime_error("hash for download does not match manifest: " + primer_url.second);
    std::ifstream f2Check(f2, std::ios::binary);
    std::vector<unsigned char> s2(picosha2::k_digest_size);
    picosha2::hash256(f2Check, s2.begin(), s2.end());
    std::string f2h;
    picosha2::bytes_to_hex_string(s2.begin(), s2.end(), f2h);
    LOG_TRACE("\tsha256 for reference:\t{}", f2h);
    if (f2h != reference_url.second)
        throw std::runtime_error("hash for download does not match manifest: " + reference_url.second);

    // pass back the downloaded scheme via the args struct
    args.schemeFile = f1;
    args.refSeqFile = f2;
}

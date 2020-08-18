#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "primerScheme.hpp"
#include "rapidcsv.h"

// some useful tags
const std::string tag_leftPrimer = "_LEFT";
const std::string tag_rightPrimer = "_RIGHT";
const std::string tag_altPrimer = "_alt";
const std::string Unmatched_Pool = "unmatched";

// Primer constructor.
artic::Primer::Primer(unsigned int start, unsigned int end, const std::string primerID, size_t poolID)
    : _start(start), _end(end), _primerID(primerID), _poolID(poolID)
{
    _numAlts = 0;

    // check the fields are valid and add the direction
    if (_primerID.empty())
        throw std::runtime_error("primer constructor received missing ID");

    // add the direction, based on the primer ID
    std::size_t left = primerID.find(tag_leftPrimer);
    std::size_t right = primerID.find(tag_rightPrimer);
    if ((left == std::string::npos) && (right == std::string::npos))
        throw std::runtime_error("invalid primer ID doesn't contain LEFT/RIGHT: " + _primerID);
    if (left != std::string::npos)
    {
        if (right != std::string::npos)
            throw std::runtime_error("invalid primer ID contains both LEFT and RIGHT: " + _primerID);
        _isForward = true;
        _baseID = _primerID.substr(0, left);
    }
    else
    {
        _isForward = false;
        _baseID = _primerID.substr(0, right);
    }

    // check that the start/end is valid
    if (_start >= _end)
        throw std::runtime_error("invalid primer start/end for primerID: " + _primerID);
}

// MergeAlt will merge a primer with an alt, yielding a primer with the maximal span.
void artic::Primer::MergeAlt(const Primer& alt)
{
    if (_isForward != alt._isForward)
        throw std::runtime_error("could not merge alt with different orientation to canonical");
    if (alt._start < _start)
        _start = alt._start;
    if (alt._end > _end)
        _end = alt._end;
    _numAlts++;
}

// GetNumAlts returns the number of alts incorporated into a primer.
unsigned int artic::Primer::GetNumAlts(void) { return _numAlts; }

// GetStart returns the primer start.
int64_t artic::Primer::GetStart(void) { return _start; }

// GetEnd returns the primer end.
int64_t artic::Primer::GetEnd(void) { return _end; }

// GetLen returns the length of the primer sequence.
unsigned int artic::Primer::GetLen(void) { return _end - _start; }

// GetID returns the primerID.
const std::string& artic::Primer::GetID(void) const { return _primerID; }

// GetBaseID returns the baseID of the primer (with no alt or LEFT/RIGHT tag).
const std::string& artic::Primer::GetBaseID(void) const { return _baseID; }

// GetPrimerPoolID returns the primer pool ID for the primer.
size_t artic::Primer::GetPrimerPoolID(void) const { return _poolID; }

// IsForward returns the primer direction (true = forward, false = reverse).
bool artic::Primer::IsForward(void) { return _isForward; }

// PrimerScheme constructor.
artic::PrimerScheme::PrimerScheme(const std::string inputFile, unsigned int schemeVersion)
    : _version(schemeVersion)
{
    // initialise the private members
    _numPrimers = 0;
    _numAlts = 0;

    // check the version provided
    if (_version < 1 || _version > 3)
        throw std::runtime_error("unrecognised primer scheme version: " + std::to_string(_version));

    // check the input file provided
    if (inputFile.empty())
        throw std::runtime_error("primer scheme input file required");

    // load the primer scheme input file and check for the expected number of rows
    std::ifstream ifile;
    ifile.open(inputFile);
    if (!ifile)
        throw std::runtime_error("primer scheme file does not exist");
    rapidcsv::Document scheme(inputFile, rapidcsv::LabelParams(-1, -1), rapidcsv::SeparatorParams('\t'));

    // TODO: this needs adjusting for different scheme versions
    if (scheme.GetColumnCount() < 5)
        throw std::runtime_error("too few columns in primer scheme file - check it's in ARTIC format");

    // iterate over the primer scheme rows
    auto rowCount = scheme.GetRowCount();
    for (unsigned int rowIterator = 0; rowIterator < rowCount; ++rowIterator)
    {
        std::vector<std::string> row = scheme.GetRow<std::string>(rowIterator);

        // check we don't have multiple references used in the scheme
        if (!_referenceID.empty())
        {
            if (row[0] != _referenceID)
                throw std::runtime_error("multiple reference sequences can't be used in primer scheme");
        }
        else
        {
            _referenceID = row[0];
        }

        // add the primer pool to the scheme and get a lookup int for the primers
        size_t poolID = 0;
        std::vector<std::string>::iterator poolItr = std::find(_primerPools.begin(), _primerPools.end(), row[4]);
        if (poolItr != _primerPools.end())
        {
            poolID = poolItr - _primerPools.begin();
        }
        else
        {
            // store the primer pool name
            _primerPools.emplace_back(row[4]);

            // add a bit vector for storing pool primer locations
            _primerSites.emplace_back(sul::dynamic_bitset());

            // get the pool ID for the new pool name
            poolID = _primerPools.size() - 1;
        }

        // try converting the primer scheme row into a primer object
        Primer* primer;
        try
        {
            primer = new Primer(std::stoi(row[1]), std::stoi(row[2]), row[3], poolID);
        }
        catch (const std::exception& e)
        {
            std::cerr << std::string("skipping row ") << rowIterator << " in scheme: " << e.what() << std::endl;
            continue;
        }

        // increment the raw primer counter
        _numPrimers++;

        // chomp off any alt tag to get the canonical primer ID
        std::string canonicalID = row[3].substr(0, row[3].find(tag_altPrimer));

        // check to see if this primer or an alt has not been seen before and then add it to the forward/reverse map
        if (primer->IsForward())
        {
            schemeMap::iterator i = _fPrimers.find(canonicalID);
            if (i == _fPrimers.end())
            {
                _fPrimers.emplace(canonicalID, primer);
                continue;
            }

            // otherwise, the primer has been seen before so it's an alt that needs merging
            i->second->MergeAlt(*primer);
            _numAlts++;
        }
        else
        {
            schemeMap::iterator j = _rPrimers.find(canonicalID);
            if (j == _rPrimers.end())
            {
                _rPrimers.emplace(canonicalID, primer);
                continue;
            }

            // otherwise, the primer has been seen before so it's an alt that needs merging
            j->second->MergeAlt(*primer);
            _numAlts++;
        }
    }
    ifile.close();

    // check that all rows were handled
    if (_numPrimers != rowCount)
        throw std::runtime_error("primer count does not equal the number of rows in the input file: " + std::to_string(_numPrimers) + " vs " + std::to_string(rowCount));

    // check the scheme
    _checkScheme();
}

// PrimerScheme destructor.
artic::PrimerScheme::~PrimerScheme(void)
{
    // properly destory the primer scheme (maps of pointers)
    for (schemeMap::iterator i = _fPrimers.begin(); i != _fPrimers.end(); ++i)
    {
        delete (i->second);
        i->second = nullptr;
    }
    for (schemeMap::iterator j = _rPrimers.begin(); j != _rPrimers.end(); ++j)
    {
        delete (j->second);
        j->second = nullptr;
    }
}

// GetVersion returns the primer scheme version.
unsigned int artic::PrimerScheme::GetVersion(void) { return _version; }

// GetReferenceID returns the reference sequence ID found in the primer scheme.
const std::string& artic::PrimerScheme::GetReferenceID(void) const { return _referenceID; }

// GetNumPrimers returns the number of primers in the primer scheme.
unsigned int artic::PrimerScheme::GetNumPrimers(void) { return _numPrimers; }

// GetNumAlts returns the number of alts in the primer scheme.
unsigned int artic::PrimerScheme::GetNumAlts(void) { return _numAlts; }

// GetNumAmplicons returns the number of primers in the primer scheme (after alts are merged).
unsigned int artic::PrimerScheme::GetNumAmplicons(void) { return _numAmplicons; }

// GetMeanAmpliconSpan returns the mean amplicon span (including primer sequence).
unsigned int artic::PrimerScheme::GetMeanAmpliconSpan(void) { return _meanAmpliconSpan; }

// GetPrimerPools returns the primer pools found in the primer scheme.
std::vector<std::string> artic::PrimerScheme::GetPrimerPools(void) { return _primerPools; }

// GetPrimerPool returns the primer pool for the provided pool ID.
const std::string& artic::PrimerScheme::GetPrimerPool(size_t poolID) const
{
    if (poolID < 0 || poolID >= _primerPools.size())
        throw std::runtime_error("poolID not found in scheme pools - " + std::to_string(poolID));
    return _primerPools.at(poolID);
}

// GetPrimerPoolID returns the primer pool ID for the provided pool name.
size_t artic::PrimerScheme::GetPrimerPoolID(const std::string& poolName) const
{
    auto itr = std::find(_primerPools.begin(), _primerPools.end(), poolName);
    if (itr != _primerPools.end())
        return itr - _primerPools.begin();
    else
        throw std::runtime_error("pool name not found in scheme - " + poolName);
}

// GetRefStart returns the first position in the reference covered by the primer scheme.
int64_t artic::PrimerScheme::GetRefStart(void) { return _refStart; }

// GetRefEnd returns the last position in the reference covered by the primer scheme.
int64_t artic::PrimerScheme::GetRefEnd(void) { return _refEnd; }

// GetNumOverlaps returns the number of reference positions covered by more than one amplicon.
unsigned int artic::PrimerScheme::GetNumOverlaps(void) { return _ampliconOverlaps.count(); }

// FindPrimers returns a primer pair with the nearest forward and reverse primer for a given segment start and end.
// Note: the primer pair may not be correctly paired, check using the IsProperlyPaired() method
artic::Amplicon artic::PrimerScheme::FindPrimers(int64_t segStart, int64_t segEnd)
{

    // get the nearest forward primer start
    std::string fPrimerID;
    std::vector<std::pair<int64_t, std::string>>::iterator fIterator;
    fIterator = std::lower_bound(_fPrimerLocations.begin(), _fPrimerLocations.end(), std::pair<int64_t, std::string>(segStart, {}));
    if (std::abs(int(fIterator->first - segStart)) <= std::abs(int(std::prev(fIterator, 1)->first - segStart)))
    {
        fPrimerID = fIterator->second;
    }
    else
    {
        fPrimerID = (std::prev(fIterator, 1))->second;
    }

    // get the nearest reverse primer end
    std::string rPrimerID;
    std::vector<std::pair<int64_t, std::string>>::iterator rIterator;
    rIterator = std::lower_bound(_rPrimerLocations.begin(), _rPrimerLocations.end(), std::pair<int64_t, std::string>(segEnd, {}));
    if (std::abs(int(rIterator->first - segEnd)) <= std::abs(int(std::prev(rIterator, 1)->first - segEnd)))
    {
        rPrimerID = rIterator->second;
    }
    else
    {
        rPrimerID = (std::prev(rIterator, 1))->second;
    }

    // lookup the primers in the scheme
    schemeMap::iterator i = _fPrimers.find(fPrimerID);
    schemeMap::iterator j = _rPrimers.find(rPrimerID);
    if ((i == _fPrimers.end()) || (j == _rPrimers.end()))
        throw std::runtime_error("primer dropped from scheme - " + fPrimerID + " & " + rPrimerID);
    return Amplicon(i->second, j->second, this);
}

// CheckAmpliconOverlap returns true if the queried position is covered by multiple primers.
bool artic::PrimerScheme::CheckAmpliconOverlap(int64_t pos)
{
    if ((_refStart > pos) || (_refEnd < pos))
        throw std::runtime_error("query position outside of primer scheme bounds");
    return _ampliconOverlaps.test(pos);
}

// CheckPrimerSite returns true if the queried position is a primer site for the given pool.
bool artic::PrimerScheme::CheckPrimerSite(int64_t pos, const std::string& poolName)
{
    if ((_refStart > pos) || (_refEnd < pos))
        throw std::runtime_error("query position outside of primer scheme bounds");
    auto poolID = GetPrimerPoolID(poolName);
    return _primerSites.at(poolID).test(pos);
}

// _checkScheme will check all forward primers have a paired reverse primer and record some primer scheme stats.
void artic::PrimerScheme::_checkScheme(void)
{
    if (_numPrimers == 0)
        throw std::runtime_error("no primers found in the provided scheme");
    if (_fPrimers.size() != _rPrimers.size())
        throw std::runtime_error("number of forward primers does not match number of reverse primers (after alt merging): " + std::to_string(_fPrimers.size()) + " vs. " + std::to_string(_rPrimers.size()));
    _numAmplicons = 0;

    // cycle through the map holding the forward primers
    int64_t spanCounter = 0;
    for (schemeMap::iterator i = _fPrimers.begin(); i != _fPrimers.end(); ++i)
    {

        // add the forward primer start position to the holder
        _fPrimerLocations.emplace_back(i->second->GetStart(), i->second->GetID());

        // find the corresponding reverse primer and add the position to the holder
        schemeMap::iterator j = _rPrimers.find(i->second->GetBaseID() + tag_rightPrimer);
        (j == _rPrimers.end()) ? throw std::runtime_error("can't find matching reverse primer for " + i->second->GetID()) : _rPrimerLocations.emplace_back(j->second->GetEnd(), j->second->GetID());

        // increment the amplicon counter and spans
        _numAmplicons++;
        spanCounter += (j->second->GetEnd() - i->second->GetStart());
    }
    _meanAmpliconSpan = spanCounter / _numAmplicons;

    // check all primers have been properly paired
    if (_numAmplicons != _fPrimers.size())
        throw std::runtime_error("number of amplicons does not match number of forward primers: " + std::to_string(_numAmplicons) + " vs " + std::to_string(_fPrimers.size()));
    if (_numAmplicons != _rPrimers.size())
        throw std::runtime_error("number of amplicons does not match number of reverse primers: " + std::to_string(_numAmplicons) + " vs " + std::to_string(_rPrimers.size()));

    // check the same number of forward and reverse primers have been collected
    if (_fPrimerLocations.size() != _rPrimerLocations.size())
        throw std::runtime_error("mismatched number of forward and reverse primer starts: " + std::to_string(_fPrimerLocations.size()) + " vs " + std::to_string(_rPrimerLocations.size()));

    // sort the start positions
    std::sort(_fPrimerLocations.begin(), _fPrimerLocations.end());
    std::sort(_rPrimerLocations.begin(), _rPrimerLocations.end());

    // update the min/max value of the scheme
    _refStart = _fPrimerLocations.front().first;
    _refEnd = _rPrimerLocations.back().first;

    // store the primer overlap regions
    _ampliconOverlaps.resize(_refEnd - _refStart, 0);
    for (unsigned int i = 0; i < _numAmplicons - 1; i++)
    {
        if (_fPrimerLocations.at(i + 1).first < _rPrimerLocations.at(i).first)
        {
            auto start = _fPrimerLocations.at(i + 1).first;
            auto len = _rPrimerLocations.at(i).first - start;
            _ampliconOverlaps.set(start, len, 1);
        }
        else
        {
            throw std::runtime_error("gap found in primer scheme: " + std::to_string(_fPrimerLocations.at(i + 1).first) + "-" + std::to_string(_rPrimerLocations.at(i).first));
        }
    }

    // store the primer sites per pool
    for (size_t poolID = 0; poolID < _primerPools.size(); ++poolID)
    {
        _primerSites.at(poolID).resize(_refEnd - _refStart, 0);
        for (auto const& primer : _fPrimers)
        {
            if (primer.second->GetPrimerPoolID() == poolID)
                _primerSites.at(poolID).set(primer.second->GetStart(), primer.second->GetLen(), 1);
        }
        for (auto const& primer : _rPrimers)
        {
            if (primer.second->GetPrimerPoolID() == poolID)
                _primerSites.at(poolID).set(primer.second->GetEnd(), primer.second->GetLen(), 1);
        }
    }
}

// Amplicon constructor.
artic::Amplicon::Amplicon(Primer* p1, Primer* p2, PrimerScheme* scheme)
    : _fPrimer(p1), _rPrimer(p2), _scheme(scheme)
{
}

// IsProperlyPaired returns true if this primer is properly paired.
bool artic::Amplicon::IsProperlyPaired(void)
{
    // paired if baseID matches, directions oppose and primer pool matches
    return (_fPrimer->GetBaseID() == _rPrimer->GetBaseID()) &&
           (_fPrimer->IsForward() != _rPrimer->IsForward()) &&
           (_fPrimer->GetPrimerPoolID() == _rPrimer->GetPrimerPoolID());
}

// GetID returns the shared ID string of the primer pair.
const std::string artic::Amplicon::GetID(void) const { return std::string(_fPrimer->GetID() + "_" + _rPrimer->GetID()); }

// GetPrimerPoolID returns the pool for the primer pair (Unmatched returned if primers not properly paired).
const std::string& artic::Amplicon::GetPrimerPool(void)
{
    if (IsProperlyPaired())
    {
        auto id = _fPrimer->GetPrimerPoolID();
        return _scheme->GetPrimerPool(id);
    }
    return Unmatched_Pool;
}

// GetMaxSpan returns the start and end of the amplicon, including the primer sequence.
std::pair<int64_t, int64_t> artic::Amplicon::GetMaxSpan(void) { return std::pair(_fPrimer->GetStart(), _rPrimer->GetEnd()); }

// GetMinSpan returns the start and end of the amplicon, excluding the primer sequence.
std::pair<int64_t, int64_t> artic::Amplicon::GetMinSpan(void) { return std::pair(_fPrimer->GetEnd(), _rPrimer->GetStart()); }

// GetForwardPrimer returns a pointer to the forward primer in the amplicon.
const artic::Primer* artic::Amplicon::GetForwardPrimer(void) { return _fPrimer; }

// GetReversePrimer returns a pointer to the reverse primer in the amplicon.
const artic::Primer* artic::Amplicon::GetReversePrimer(void) { return _rPrimer; }
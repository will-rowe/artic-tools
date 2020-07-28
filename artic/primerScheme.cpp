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
artic::Primer::Primer(const std::string chrom, unsigned int start, unsigned int end, const std::string primerID, const std::string poolName)
    : _chrom(chrom), _start(start), _end(end), _primerID(primerID), _poolName(poolName)
{
    // initialise some members
    _numAlts = 0;
    _direction = 0;

    // check the fields are valid and add the direction
    if (_chrom.empty() || _primerID.empty() || _poolName.empty())
    {
        throw std::runtime_error("primer constructor received missing field");
    }

    // add the direction, based on the primer ID
    std::size_t left = primerID.find(tag_leftPrimer);
    std::size_t right = primerID.find(tag_rightPrimer);
    if ((left == std::string::npos) && (right == std::string::npos))
    {
        throw std::runtime_error("invalid primer ID doesn't contain LEFT/RIGHT: " + _primerID);
    }
    if (left != std::string::npos)
    {
        if (right != std::string::npos)
        {
            throw std::runtime_error("invalid primer ID contains both LEFT and RIGHT: " + _primerID);
        }
        _direction = 1;
        _baseID = _primerID.substr(0, left);
    }
    else
    {
        _direction = -1;
        _baseID = _primerID.substr(0, right);
    }

    // check that the start/end is valid
    if (_start >= _end)
    {
        throw std::runtime_error("invalid primer start/end for primerID: " + _primerID);
    }
}

// MergeAlt will merge a primer with an alt, yielding a primer with the maximal span.
void artic::Primer::MergeAlt(const Primer& alt)
{
    if (_direction != alt._direction)
    {
        throw std::runtime_error("could not merge alt with different orientation to canonical");
    }
    if (alt._start < _start)
    {
        _start = alt._start;
    }
    if (alt._end > _end)
    {
        _end = alt._end;
    }
    _numAlts++;
}

// GetNumAlts returns the number of alts incorporated into a primer.
unsigned int artic::Primer::GetNumAlts(void) { return _numAlts; }

// GetStart returns the primer start.
unsigned int artic::Primer::GetStart(void) { return _start; }

// GetEnd returns the primer end.
unsigned int artic::Primer::GetEnd(void) { return _end; }

// GetID returns the primerID.
const std::string& artic::Primer::GetID(void) const { return _primerID; }

// GetDirection returns the primer direction.
signed int artic::Primer::GetDirection(void) { return _direction; }

// GetBaseID returns the baseID of the primer (with no alt or LEFT/RIGHT tag).
const std::string& artic::Primer::GetBaseID(void) { return _baseID; }

// GetPrimerPool returns the primer pool for the primer.
const std::string& artic::Primer::GetPrimerPool(void) { return _poolName; }

// Amplicon constructor.
artic::Amplicon::Amplicon(Primer* p1, Primer* p2)
    : _fPrimer(p1), _rPrimer(p2)
{
}

// IsProperlyPaired returns true if this primer is properly paired.
bool artic::Amplicon::IsProperlyPaired(void)
{
    // paired if baseID matches, directions oppose and primer pool matches
    return (_fPrimer->GetBaseID() == _rPrimer->GetBaseID()) &&
           (_fPrimer->GetDirection() != _rPrimer->GetDirection()) &&
           (_fPrimer->GetPrimerPool() == _rPrimer->GetPrimerPool());
}

// GetID returns the shared ID string of the primer pair.
const std::string artic::Amplicon::GetID(void)
{
    return std::string(_fPrimer->GetID() + "_" + _rPrimer->GetID());
}

// GetPrimerPool returns the pool for the primer pair (Unmatched_Pool if not properly paired).
const std::string& artic::Amplicon::GetPrimerPool(void)
{
    if (IsProperlyPaired())
    {
        return _fPrimer->GetPrimerPool();
    }
    return Unmatched_Pool;
}

// GetMaxSpan returns the start and end of the amplicon, including the primer sequence.
std::pair<unsigned int, unsigned int> artic::Amplicon::GetMaxSpan(void)
{
    std::pair<unsigned int, unsigned int> span;
    span.first = _fPrimer->GetStart();
    span.second = _rPrimer->GetEnd();
    return span;
}

// GetMinSpan returns the start and end of the amplicon, excluding the primer sequence.
std::pair<unsigned int, unsigned int> artic::Amplicon::GetMinSpan(void)
{
    std::pair<unsigned int, unsigned int> span;
    span.first = _fPrimer->GetEnd();
    span.second = _rPrimer->GetStart();
    return span;
}

// GetForwardPrimer returns a pointer to the forward primer in the amplicon.
const artic::Primer* artic::Amplicon::GetForwardPrimer(void)
{
    return _fPrimer;
}

// GetReversePrimer returns a pointer to the reverse primer in the amplicon.
const artic::Primer* artic::Amplicon::GetReversePrimer(void)
{
    return _rPrimer;
}

// PrimerScheme constructor.
artic::PrimerScheme::PrimerScheme(const std::string inputFile, unsigned int schemeVersion)
    : _version(schemeVersion)
{
    // initialise the private members
    _numPrimers = 0;
    _numAlts = 0;
    _primerPools.emplace_back(Unmatched_Pool);

    // check the version provided
    if (_version < 1 || _version > 3)
    {
        throw std::runtime_error("unrecognised primer scheme version: " + std::to_string(_version));
    }

    // check the input file provided
    if (inputFile.empty())
    {
        throw std::runtime_error("primer scheme input file required");
    }

    // load the primer scheme input file and check for the expected number of rows
    std::ifstream ifile;
    ifile.open(inputFile);
    if (!ifile)
    {
        throw std::runtime_error("primer scheme file does not exist");
    }
    rapidcsv::Document scheme(inputFile, rapidcsv::LabelParams(-1, -1), rapidcsv::SeparatorParams('\t'));
    if (scheme.GetColumnCount() < 5)
    {
        throw std::runtime_error("too few columns in primer scheme file - check it's in ARTIC format");
    }

    // iterate over the primer scheme rows
    auto rowCount = scheme.GetRowCount();
    for (unsigned int rowIterator = 0; rowIterator < rowCount; ++rowIterator)
    {

        // try converting the primer scheme row into a primer object
        std::vector<std::string> row = scheme.GetRow<std::string>(rowIterator);
        Primer* primer;
        try
        {
            primer = new Primer(row[0], std::stoi(row[1]), std::stoi(row[2]), row[3], row[4]);
        }
        catch (const std::exception& e)
        {
            std::cerr << e.what() << std::string(" - skipping row: ") << rowIterator << std::endl;
            continue;
        }

        // increment the raw primer counter
        _numPrimers++;

        // add the primer pool
        std::vector<std::string>::iterator it = std::find(_primerPools.begin(), _primerPools.end(), primer->GetPrimerPool());
        if (it == _primerPools.end())
        {
            _primerPools.emplace_back(primer->GetPrimerPool());
        }

        // chomp off any alt tag to get the canonical primer ID
        std::string canonicalID = row[3].substr(0, row[3].find(tag_altPrimer));

        // check to see if this primer or an alt has not been seen before and then add it to the forward/reverse map
        if (primer->GetDirection() == 1)
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
    {
        throw std::runtime_error("primer count does not equal the number of rows in the input file: " + std::to_string(_numPrimers) + " vs " + std::to_string(rowCount));
    }

    // check the scheme
    CheckScheme();
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

// CheckScheme will check all forward primers have a paired reverse primer and record some primer scheme stats.
void artic::PrimerScheme::CheckScheme(void)
{
    if (_numPrimers == 0)
    {
        throw std::runtime_error("no primers found in the provided scheme");
    }
    if (_fPrimers.size() != _rPrimers.size())
    {
        throw std::runtime_error("number of forward primers does not match number of reverse primers (after alt merging)");
    }
    _numAmplicons = 0;

    // cycle through the map holding the forward primers
    for (schemeMap::iterator i = _fPrimers.begin(); i != _fPrimers.end(); ++i)
    {
        // add the forward primer start position to the holder
        _fPrimerLocations.emplace_back(i->second->GetStart(), i->second->GetID());

        // find the corresponding reverse primer
        schemeMap::iterator j = _rPrimers.find(i->second->GetBaseID() + tag_rightPrimer);
        if (j == _rPrimers.end())
        {
            throw std::runtime_error("can't find matching reverse primer for " + i->second->GetID());
        }

        // add the reverse primer end position to the holder
        _rPrimerLocations.emplace_back(j->second->GetEnd(), j->second->GetID());
        _numAmplicons++;
    }

    // check all primers have been properly paired
    if (_numAmplicons != _fPrimers.size())
    {
        throw std::runtime_error("number of amplicons does not match number of forward primers: " + std::to_string(_numAmplicons) + " vs " + std::to_string(_fPrimers.size()));
    }

    // sort the start positions
    std::sort(_fPrimerLocations.begin(), _fPrimerLocations.end());
    std::sort(_rPrimerLocations.begin(), _rPrimerLocations.end());

    // check the same number of forward and reverse primers have been collected
    if (_fPrimerLocations.size() != _rPrimerLocations.size())
    {
        throw std::runtime_error("primer scheme contains mismatched number of forward and reverse primers (post merging alts): " + std::to_string(_fPrimerLocations.size()) +
                                 " vs " +
                                 std::to_string(_rPrimerLocations.size()));
    }

    // update the min/max value of the scheme
    _minStart = _fPrimerLocations.front().first;
    _maxEnd = _rPrimerLocations.back().first;
}

// GetVersion returns the primer scheme version.
unsigned int artic::PrimerScheme::GetVersion(void) { return _version; }

// GetNumPrimers returns the number of primers in the primer scheme.
unsigned int artic::PrimerScheme::GetNumPrimers(void) { return _numPrimers; }

// GetNumAlts returns the number of alts in the primer scheme.
unsigned int artic::PrimerScheme::GetNumAlts(void) { return _numAlts; }

// GetNumAmplicons returns the number of primers in the primer scheme (after alts are merged).
unsigned int artic::PrimerScheme::GetNumAmplicons(void) { return _numAmplicons; }

// GetPrimerPools returns the primer pool found in the primer scheme
std::vector<std::string> artic::PrimerScheme::GetPrimerPools(void) { return _primerPools; }

// FindPrimers returns a primer pair with the nearest forward and reverse primer for a given segment start and end.
// Note: the primer pair may not be correctly paired, check using the IsProperlyPaired() method
artic::Amplicon artic::PrimerScheme::FindPrimers(unsigned int segStart, unsigned int segEnd)
{
    //if ((_minStart > segStart) || (_maxEnd < segEnd))
    //{
    //    throw std::runtime_error("alignment is outside of available primer scheme bounds");
    //}

    // get the nearest forward primer start
    std::string fPrimerID;
    std::vector<std::pair<unsigned int, std::string>>::iterator fIterator;
    fIterator = std::lower_bound(_fPrimerLocations.begin(), _fPrimerLocations.end(), std::pair<unsigned int, std::string>(segStart, {}));
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
    std::vector<std::pair<unsigned int, std::string>>::iterator rIterator;
    rIterator = std::lower_bound(_rPrimerLocations.begin(), _rPrimerLocations.end(), std::pair<unsigned int, std::string>(segEnd, {}));
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
    {
        throw std::runtime_error("primer dropped from scheme - " + fPrimerID + " & " + rPrimerID);
    }
    return Amplicon(i->second, j->second);
}

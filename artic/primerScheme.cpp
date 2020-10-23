#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <rapidcsv.h>
#include <stdexcept>
#include <string>

#include "primerScheme.hpp"

// ARTIC scheme tags
const std::string LEFT_PRIMER_TAG = "_LEFT";
const std::string RIGHT_PRIMER_TAG = "_RIGHT";
const std::string ALT_PRIMER_TAG = "_alt";
const std::string NO_POOL = "unmatched";

// Primer constructor.
artic::Primer::Primer(unsigned int start, unsigned int end, const std::string primerID, size_t poolID)
    : _start(start), _end(end), _primerID(primerID), _poolID(poolID)
{
    _numAlts = 0;

    // check the fields are valid and add the direction
    if (_primerID.empty())
        throw std::runtime_error("primer constructor received missing ID");
    if (_start >= _end)
        throw std::runtime_error("invalid primer start/end for primerID: " + _primerID);

    // add the direction, based on the primer ID
    std::size_t left = primerID.find(LEFT_PRIMER_TAG);
    std::size_t right = primerID.find(RIGHT_PRIMER_TAG);
    if ((left == std::string::npos) && (right == std::string::npos))
        throw std::runtime_error("invalid primer ID doesn't contain LEFT/RIGHT: " + _primerID);
    if (left != std::string::npos)
    {
        if (right != std::string::npos)
            throw std::runtime_error("invalid primer ID contains both LEFT and RIGHT: " + _primerID);
        _isForward = true;
        _baseIDit = left;
        return;
    }
    _isForward = false;
    _baseIDit = right;
}

// MergeAlt will merge a primer with an alt, yielding a primer with the maximal span.
void artic::Primer::MergeAlt(const Primer& alt)
{
    if (_isForward != alt._isForward)
        throw std::runtime_error("could not merge alt with different orientation to canonical");
    if (_poolID != alt.GetPrimerPoolID())
        throw std::runtime_error("could not merge alt from different pool to canonical");
    if (alt._start < _start)
        _start = alt._start;
    if (alt._end > _end)
        _end = alt._end;
    _numAlts++;
}

// GetNumAlts returns the number of alts incorporated into a primer.
unsigned int artic::Primer::GetNumAlts(void) const { return _numAlts; }

// GetStart returns the primer start.
int64_t artic::Primer::GetStart(void) const { return _start; }

// GetEnd returns the primer end.
int64_t artic::Primer::GetEnd(void) const { return _end; }

// GetLen returns the length of the primer sequence.
unsigned int artic::Primer::GetLen(void) const { return _end - _start; }

// GetID returns the primerID.
const std::string& artic::Primer::GetID(void) const { return _primerID; }

// GetBaseID returns the baseID of the primer (with no alt or LEFT/RIGHT tag).
std::string artic::Primer::GetBaseID(void) const { return _primerID.substr(0, _baseIDit); }

// GetPrimerPoolID returns the primer pool ID for the primer.
size_t artic::Primer::GetPrimerPoolID(void) const { return _poolID; }

// IsForward returns the primer direction (true = forward, false = reverse).
bool artic::Primer::IsForward(void) { return _isForward; }

// GetSeq returns the primer sequence from a reference.
const std::string artic::Primer::GetSeq(faidx_t* reference, const std::string& refID) const
{
    if (!reference)
        throw std::runtime_error("no reference fasta provided");
    int len;
    char* seq = faidx_fetch_seq(
        reference,
        refID.c_str(),
        _start,
        _end - 1,
        &len);
    if (!seq)
        throw std::runtime_error("cannot fetch the reference sequence");
    if (len != int(GetLen()))
        throw std::runtime_error("did not fetch correct number of primer bases (got " + std::to_string(len) + ")");
    return seq;
}

// PrimerScheme constructor.
artic::PrimerScheme::PrimerScheme(const std::string& inputFile)
    : _filename(inputFile)
{
    _loadScheme(_filename);
    _validateScheme();
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

// GetFileName returns the filename that the primer scheme was loaded from.
const std::string& artic::PrimerScheme::GetFileName(void) const { return _filename; }

// GetReferenceName returns the reference sequence ID found in the primer scheme.
const std::string& artic::PrimerScheme::GetReferenceName(void) const { return _referenceID; }

// GetNumPrimers returns the number of primers in the primer scheme.
unsigned int artic::PrimerScheme::GetNumPrimers(void) { return _numPrimers; }

// GetNumAlts returns the number of alts in the primer scheme.
unsigned int artic::PrimerScheme::GetNumAlts(void) { return _numAlts; }

// GetNumAmplicons returns the number of primers in the primer scheme (after alts are merged).
unsigned int artic::PrimerScheme::GetNumAmplicons(void) { return _numAmplicons; }

// GetMeanAmpliconSpan returns the mean amplicon span (including primer sequence).
unsigned int artic::PrimerScheme::GetMeanAmpliconSpan(void) { return _meanAmpliconSpan; }

// GetPrimerPools returns the primer pools found in the primer scheme.
std::vector<std::string> artic::PrimerScheme::GetPrimerPools(void) { return std::vector<std::string>(_primerPools.begin() + 1, _primerPools.end()); }

// GetPrimerPool returns the primer pool for the provided pool ID.
const std::string& artic::PrimerScheme::GetPrimerPool(size_t poolID) const
{
    if (poolID >= _primerPools.size())
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

// GetExpAmplicons returns a vector to the amplicons the scheme expects to produce.
const std::vector<artic::Amplicon>& artic::PrimerScheme::GetExpAmplicons(void)
{
    std::sort(_expAmplicons.begin(), _expAmplicons.end(), [](auto& lhs, auto& rhs) {
        return lhs.GetForwardPrimer()->GetEnd() < rhs.GetForwardPrimer()->GetEnd();
    });
    return _expAmplicons;
}

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
    return Amplicon(i->second, j->second);
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
    return _primerSites.test(pos + (_refEnd * poolID));
}

// _loadScheme will load an input file and create the primer objects.
void artic::PrimerScheme::_loadScheme(const std::string& filename)
{
    _numPrimers = 0;
    _numAlts = 0;
    _primerPools.emplace_back(NO_POOL);

    // check the input file is provided
    if (filename.empty())
        throw std::runtime_error("primer scheme input file required");

    // load the primer scheme input file and check for the expected number of rows
    std::ifstream ifile;
    ifile.open(filename);
    if (!ifile)
        throw std::runtime_error("primer scheme file does not exist");
    rapidcsv::Document scheme(filename, rapidcsv::LabelParams(-1, -1), rapidcsv::SeparatorParams('\t'));

    // there needs to be at least 5 rows to check or it won't be a scheme we can use
    if (scheme.GetColumnCount() < 5)
        throw std::runtime_error("less than 5 columns in the primer scheme file - check it's in ARTIC format");

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
            poolID = _primerPools.size() - 1;
        }

        // try converting the primer scheme row into a primer object
        Primer* primer = 0;
        try
        {
            primer = new Primer(std::stoi(row[1]), std::stoi(row[2]), row[3], poolID);
        }
        catch (const std::exception& e)
        {
            std::cerr << "skipping row " << rowIterator << " in scheme - " << e.what() << std::endl;
            continue;
        }
        if (!primer)
            std::cerr << "failed to produce primer from row " << rowIterator << " in primer scheme" << std::endl;

        // increment the raw primer counter
        _numPrimers++;

        // chomp off any alt tag to get the canonical primer ID
        std::string canonicalID = row[3].substr(0, row[3].find(ALT_PRIMER_TAG));

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
        throw std::runtime_error("primer count does not equal the number of rows in the input file - " + std::to_string(_numPrimers) + " vs " + std::to_string(rowCount));
}

// _validateScheme will check all forward primers have a paired reverse primer and record some primer scheme stats.
void artic::PrimerScheme::_validateScheme(void)
{
    if (_numPrimers == 0)
        throw std::runtime_error("no primers found in the provided scheme");
    if (_fPrimers.size() != _rPrimers.size())
        throw std::runtime_error("number of forward primers does not match number of reverse primers (after alt merging) - " + std::to_string(_fPrimers.size()) + " vs. " + std::to_string(_rPrimers.size()));
    _numAmplicons = 0;

    // cycle through the map holding the forward primers
    int64_t spanCounter = 0;
    for (schemeMap::iterator i = _fPrimers.begin(); i != _fPrimers.end(); ++i)
    {

        // add the forward primer start position to the holder
        _fPrimerLocations.emplace_back(i->second->GetStart(), i->second->GetID());

        // find the corresponding reverse primer and add the position to the holder
        schemeMap::iterator j = _rPrimers.find(i->second->GetBaseID() + RIGHT_PRIMER_TAG);
        (j == _rPrimers.end()) ? throw std::runtime_error("can't find matching reverse primer for " + i->second->GetID()) : _rPrimerLocations.emplace_back(j->second->GetEnd(), j->second->GetID());

        // increment the amplicon counter and spans
        _numAmplicons++;
        spanCounter += (j->second->GetEnd() - i->second->GetStart());

        // create an amplicon and add it to the scheme holder
        _expAmplicons.emplace_back(Amplicon(i->second, j->second));
    }
    _meanAmpliconSpan = spanCounter / _numAmplicons;

    // check all primers have been properly paired
    if (_numAmplicons != _fPrimers.size())
        throw std::runtime_error("number of amplicons does not match number of forward primers - " + std::to_string(_numAmplicons) + " vs " + std::to_string(_fPrimers.size()));
    if (_numAmplicons != _rPrimers.size())
        throw std::runtime_error("number of amplicons does not match number of reverse primers - " + std::to_string(_numAmplicons) + " vs " + std::to_string(_rPrimers.size()));

    // check the same number of forward and reverse primers have been collected
    if (_fPrimerLocations.size() != _rPrimerLocations.size())
        throw std::runtime_error("mismatched number of forward and reverse primer starts - " + std::to_string(_fPrimerLocations.size()) + " vs " + std::to_string(_rPrimerLocations.size()));

    // sort the start positions
    std::sort(_fPrimerLocations.begin(), _fPrimerLocations.end());
    std::sort(_rPrimerLocations.begin(), _rPrimerLocations.end());

    // update the min/max value of the scheme
    _refStart = _fPrimerLocations.front().first;
    _refEnd = _rPrimerLocations.back().first;

    // store the primer overlap regions
    _ampliconOverlaps.resize(_refEnd, 0);
    for (unsigned int i = 0; i < _numAmplicons - 1; i++)
    {
        if (_fPrimerLocations.at(i + 1).first < _rPrimerLocations.at(i).first)
        {
            for (auto bitSetter = _fPrimerLocations.at(i + 1).first; bitSetter < _rPrimerLocations.at(i).first; bitSetter++)
                _ampliconOverlaps[bitSetter] = 1;
        }
        else
        {
            throw std::runtime_error("gap found in primer scheme - " + std::to_string(_fPrimerLocations.at(i + 1).first) + "-" + std::to_string(_rPrimerLocations.at(i).first));
        }
    }

    // store the primer sites per pool
    _primerSites.resize(_refEnd * _primerPools.size(), 0);
    for (size_t poolID = 0; poolID < _primerPools.size(); ++poolID)
    {
        for (auto const& primer : _fPrimers)
        {
            if (primer.second->GetPrimerPoolID() == poolID)
            {

                for (auto bitSetter = (primer.second->GetStart() + (_refEnd * poolID)); bitSetter < (primer.second->GetEnd() + (_refEnd * poolID)); bitSetter++)
                    _primerSites[bitSetter] = 1;
            }
        }
        for (auto const& primer : _rPrimers)
        {
            if (primer.second->GetPrimerPoolID() == poolID)
            {
                for (auto bitSetter = (primer.second->GetEnd() + (_refEnd * poolID)); bitSetter < (primer.second->GetStart() + (_refEnd * poolID)); bitSetter++)
                    _primerSites[bitSetter] = 1;
            }
        }
    }
}

// Amplicon constructor.
artic::Amplicon::Amplicon(Primer* p1, Primer* p2)
    : _fPrimer(p1), _rPrimer(p2)
{
    // ensure p1 is forward and p2 is reverse
    if (_fPrimer->IsForward() == _rPrimer->IsForward())
        throw std::runtime_error("cannot create amplicon from primers with the same directionality");
    if (!_fPrimer->IsForward())
    {
        _fPrimer = p2;
        _rPrimer = p1;
    }

    // p1 must come before p2
    if (_fPrimer->GetEnd() >= _rPrimer->GetStart())
        throw std::runtime_error("cannnot create amplicon from outward facing primers");

    // set if properly paired
    _isProperlyPaired = (_fPrimer->GetBaseID() == _rPrimer->GetBaseID()) &&
                        (_fPrimer->GetPrimerPoolID() == _rPrimer->GetPrimerPoolID());
}

// IsProperlyPaired returns true if this primer is properly paired.
bool artic::Amplicon::IsProperlyPaired(void) { return _isProperlyPaired; }

// GetID returns the shared ID string of the primer pair.
const std::string artic::Amplicon::GetID(void) const { return std::string(_fPrimer->GetID() + "_" + _rPrimer->GetID()); }

// GetPrimerPoolID returns the pool ID for the primer pair (0 returned if primers not properly paired).
std::size_t artic::Amplicon::GetPrimerPoolID(void)
{
    if (!_isProperlyPaired)
    {
        return 0;
    }
    return _fPrimer->GetPrimerPoolID();
}

// GetMaxSpan returns the start and end of the amplicon, including the primer sequence.
std::pair<int64_t, int64_t> artic::Amplicon::GetMaxSpan(void) { return std::pair(_fPrimer->GetStart(), _rPrimer->GetEnd()); }

// GetMinSpan returns the start and end of the amplicon, excluding the primer sequence.
std::pair<int64_t, int64_t> artic::Amplicon::GetMinSpan(void) { return std::pair(_fPrimer->GetEnd(), _rPrimer->GetStart()); }

// GetForwardPrimer returns a pointer to the forward primer in the amplicon.
const artic::Primer* artic::Amplicon::GetForwardPrimer(void) { return _fPrimer; }

// GetReversePrimer returns a pointer to the reverse primer in the amplicon.
const artic::Primer* artic::Amplicon::GetReversePrimer(void) { return _rPrimer; }

// AddKmers adds the k-mers from a sequence to the amplicon.
void artic::Amplicon::AddKmers(const char* seq, uint32_t seqLen, uint32_t kSize)
{
    artic::getEncodedKmers(seq, seqLen, kSize, _kmers);
    return;
}
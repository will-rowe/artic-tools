#ifndef PRIMERSCHEME_H
#define PRIMERSCHEME_H

#include <boost/dynamic_bitset.hpp>
#include <htslib/faidx.h>
#include <string>
#include <unordered_map>
#include <vector>

namespace artic
{
    class Primer;
    class PrimerScheme;
    class Amplicon;
    typedef std::unordered_map<std::string, Primer*> schemeMap;

    //******************************************************************************
    // Primer class holds primer information.
    //
    // NOTES:
    // * primer direction assumes that 'LEFT' or 'RIGHT' is included in the primerID
    // * primer alt assumes that '_alt' is included in the primerID once only
    // * any primers with duplicate names will be merged
    //******************************************************************************
    class Primer
    {
    public:
        // Primer constructor.
        Primer(unsigned int start, unsigned int end, const std::string primerID, size_t poolID);

        // MergeAlt will merge a primer with an alt, yielding a primer with the maximal span.
        void MergeAlt(const Primer& alt);

        // GetNumAlts returns the number of alts incorporated into this primer.
        unsigned int GetNumAlts(void);

        // GetStart returns the primer start.
        int64_t GetStart(void) const;

        // GetEnd returns the primer end.
        int64_t GetEnd(void) const;

        // GetLen returns the length of the primer sequence.
        unsigned int GetLen(void) const;

        // GetID returns the primerID.
        const std::string& GetID(void) const;

        // GetBaseID returns the baseID of the primer (with _LEFT/_RIGHT stripped).
        const std::string& GetBaseID(void) const;

        // GetPrimerPoolID returns the pool ID for the primer pair.
        size_t GetPrimerPoolID(void) const;

        // IsForward returns the primer direction (true = forward, false = reverse).
        bool IsForward(void);

        // GetSeq returns the primer sequence from a reference.
        const std::string GetSeq(faidx_t* reference, const std::string& refID) const;

    private:
        int64_t _start;
        int64_t _end;
        std::string _primerID;
        size_t _poolID;
        bool _isForward;
        unsigned int _numAlts;
        std::string _baseID;
    };

    //******************************************************************************
    // PrimerScheme class handles the ARTIC style primer schemes.
    //
    //  NOTES:
    //  * only V3 tested so far
    //******************************************************************************
    class PrimerScheme
    {
    public:
        // PrimerScheme constructor and destructor.
        PrimerScheme(const std::string inputFile, unsigned int schemeVersion);
        ~PrimerScheme(void);

        // GetVersion returns the ARTIC primer scheme version (1/2/3).
        unsigned int GetVersion(void);

        // GetReferenceName returns the reference sequence ID found in the primer scheme.
        const std::string& GetReferenceName(void) const;

        // GetNumPrimers returns the total number of primers in the scheme.
        unsigned int GetNumPrimers(void);

        // GetNumAlts returns the number of alts in the primer scheme.
        unsigned int GetNumAlts(void);

        // GetNumAmplicons returns the number of amplicons in the primer scheme (after alts merged and primer pairs matched).
        unsigned int GetNumAmplicons(void);

        // GetMeanAmpliconSpan returns the mean amplicon span (including primer sequence).
        unsigned int GetMeanAmpliconSpan(void);

        // GetPrimerPools returns all the primer pools found in the primer scheme.
        std::vector<std::string> GetPrimerPools(void);

        // GetPrimerPool returns the primer pool for the provided pool ID.
        const std::string& GetPrimerPool(size_t poolID) const;

        // GetPrimerPoolID returns the primer pool ID for the provided pool name.
        size_t GetPrimerPoolID(const std::string& poolName) const;

        // GetRefStart returns the first position in the reference covered by the primer scheme.
        int64_t GetRefStart(void);

        // GetRefEnd returns the last position in the reference covered by the primer scheme.
        int64_t GetRefEnd(void);

        // GetNumOverlaps returns the number of reference positions covered by more than one amplicon.
        unsigned int GetNumOverlaps(void);

        // FindPrimers returns pointers to the nearest forward and reverse primer, given an alignment segment's start and end position.
        Amplicon FindPrimers(int64_t segStart, int64_t segEnd);

        // CheckAmpliconOverlap returns true if the queried position is covered by multiple amplicons (incl. primer sequence).
        bool CheckAmpliconOverlap(int64_t pos);

        // CheckPrimerSite returns true if the queried position is a primer site for the given pool.
        bool CheckPrimerSite(int64_t pos, const std::string& poolName);

    private:
        void _checkScheme(void); // _checkScheme will check all forward primers have a paired reverse primer and record some primer scheme stats

        unsigned int _version;                                          // the primer scheme version (based on the ARTIC versioning)
        std::string _referenceID;                                       // the ID of the reference sequence covered by the primer scheme
        unsigned int _numPrimers;                                       // the total number of primers in the scheme
        unsigned int _numAlts;                                          // the number of alts that were merged when the scheme was read
        unsigned int _numAmplicons;                                     // the number of amplicons in the scheme
        unsigned int _meanAmpliconSpan;                                 // the mean amplicon span
        int64_t _refStart;                                              // the first position in the reference covered by the primer scheme
        int64_t _refEnd;                                                // the last position in the reference covered by the primer scheme
        std::vector<std::string> _primerPools;                          // the primer pool IDs found in the primer scheme
        schemeMap _fPrimers;                                            // the forward primers for the scheme
        schemeMap _rPrimers;                                            // the reverse primers for the scheme
        std::vector<std::pair<int64_t, std::string>> _fPrimerLocations; // the start position and primerID of each forward primer in the scheme
        std::vector<std::pair<int64_t, std::string>> _rPrimerLocations; // the end position and primerID of each reverse primer in the scheme
        boost::dynamic_bitset<> _ampliconOverlaps;                      // bit vector encoding all the overlap positions in the scheme
        boost::dynamic_bitset<> _primerSites;                           // primer sites, stored in a bit vector and offset by primer pool ID
        //std::unordered_map<std::string, Amplicon> _amplicons;         // the expected amplicons produced by the scheme
    };

    //******************************************************************************
    // Amplicon is a container for two primers.
    //******************************************************************************
    class Amplicon
    {
    public:
        // Amplicon constructor.
        Amplicon(Primer* p1, Primer* p2, PrimerScheme* scheme);

        // IsProperlyPaired returns true if the amplicon primers are properly paired.
        bool IsProperlyPaired(void);

        // GetID returns the ID string of the primer pair (combines primer IDs).
        const std::string GetID(void) const;

        // GetPrimerPool returns the pool for the primer pair (Unmatched returned if not properly paired).
        const std::string& GetPrimerPool(void);

        // GetMaxSpan returns the start and end of the amplicon, including the primer sequence.
        std::pair<int64_t, int64_t> GetMaxSpan(void);

        // GetMinSpan returns the start and end of the amplicon, excluding the primer sequence.
        std::pair<int64_t, int64_t> GetMinSpan(void);

        // GetForwardPrimer returns a pointer to the forward primer in the amplicon.
        const Primer* GetForwardPrimer(void);

        // GetReversePrimer returns a pointer to the reverse primer in the amplicon.
        const Primer* GetReversePrimer(void);

    private:
        Primer* _fPrimer;      // pointer to the forward primer object
        Primer* _rPrimer;      // pointer to the reverse primer object
        PrimerScheme* _scheme; // pointer to the primer scheme
    };

} // namespace artic

#endif
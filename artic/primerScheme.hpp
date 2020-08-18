#ifndef PRIMERSCHEME_H
#define PRIMERSCHEME_H

#include <string>
#include <unordered_map>
#include <vector>

#include <sul/dynamic_bitset.hpp>

namespace artic
{

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
        Primer(unsigned int start, unsigned int end, const std::string primerID, const std::string poolName);

        // MergeAlt will merge a primer with an alt, yielding a primer with the maximal span.
        void MergeAlt(const Primer& alt);

        // GetNumAlts returns the number of alts incorporated into this primer.
        unsigned int GetNumAlts(void);

        // GetStart returns the primer start.
        int64_t GetStart(void);

        // GetEnd returns the primer end.
        int64_t GetEnd(void);

        // GetLen returns the length of the primer sequence.
        unsigned int GetLen(void);

        // GetID returns the primerID.
        const std::string& GetID(void) const;

        // GetBaseID returns the baseID of the primer (with _LEFT/_RIGHT stripped).
        const std::string& GetBaseID(void) const;

        // GetPrimerPool returns the pool for the primer pair.
        const std::string& GetPrimerPool(void) const;

        // IsForward returns the primer direction (true = forward, false = reverse).
        bool IsForward(void);

    private:
        int64_t _start;
        int64_t _end;
        std::string _primerID;
        std::string _poolName;
        bool _isForward;
        unsigned int _numAlts;
        std::string _baseID;
    };

    // schemeMap is a typedefd map for holding primerID to primer object pointer
    typedef std::unordered_map<std::string, Primer*> schemeMap;

    //******************************************************************************
    // Amplicon is a container for two primers.
    //******************************************************************************
    class Amplicon
    {
    public:
        // Amplicon constructor.
        Amplicon(Primer* p1, Primer* p2);

        // IsProperlyPaired returns true if the amplicon primers are properly paired.
        bool IsProperlyPaired(void);

        // GetID returns the ID string of the primer pair (combines primer IDs).
        const std::string GetID(void);

        // GetPrimerPool returns the pool for the primer pair (Unmatched_Pool if not properly paired).
        const std::string& GetPrimerPool(void);

        // GetMaxSpan returns the start and end of the amplicon, including the primer sequence.
        std::pair<int64_t, int64_t> GetMaxSpan(void);

        // GetMinSpan returns the start and end of the amplicon, excluding the primer sequence.
        std::pair<int64_t, int64_t> GetMinSpan(void);

        // GetForwardPrimer returns a reference to the forward primer in the amplicon.
        const Primer* GetForwardPrimer(void);

        // GetReversePrimer returns a reference to the reverse primer in the amplicon.
        const Primer* GetReversePrimer(void);

    private:
        Primer* _fPrimer; // pointer to the forward primer object
        Primer* _rPrimer; // pointer to the reverse primer object
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

        // GetReferenceID returns the reference sequence ID found in the primer scheme.
        const std::string& GetReferenceID(void) const;

        // GetNumPrimers returns the total number of primers in the scheme.
        unsigned int GetNumPrimers(void);

        // GetNumAlts returns the number of alts in the primer scheme.
        unsigned int GetNumAlts(void);

        // GetNumAmplicons returns the number of amplicons in the primer scheme (after alts merged and primer pairs matched).
        unsigned int GetNumAmplicons(void);

        // GetMeanAmpliconSpan returns the mean amplicon span (including primer sequence).
        unsigned int GetMeanAmpliconSpan(void);

        // GetPrimerPools returns the primer pools found in the primer scheme.
        std::vector<std::string> GetPrimerPools(void);

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
        bool CheckPrimerSite(int64_t pos, std::string pool);

    private:
        void _checkScheme(void); // _checkScheme will check all forward primers have a paired reverse primer and record some primer scheme stats

        unsigned int _version;                                               // the primer scheme version (based on the ARTIC versioning)
        std::string _referenceID;                                            // the ID of the reference sequence covered by the primer scheme
        unsigned int _numPrimers;                                            // the total number of primers in the scheme
        unsigned int _numAlts;                                               // the number of alts that were merged when the scheme was read
        unsigned int _numAmplicons;                                          // the number of amplicons in the scheme
        unsigned int _meanAmpliconSpan;                                      // the mean amplicon span
        int64_t _refStart;                                                   // the first position in the reference covered by the primer scheme
        int64_t _refEnd;                                                     // the last position in the reference covered by the primer scheme
        std::vector<std::string> _primerPools;                               // the primer pool IDs found in the primer scheme
        schemeMap _fPrimers;                                                 // the forward primers for the scheme
        schemeMap _rPrimers;                                                 // the reverse primers for the scheme
        std::vector<std::pair<int64_t, std::string>> _fPrimerLocations;      // the start position and primerID of each forward primer in the scheme
        std::vector<std::pair<int64_t, std::string>> _rPrimerLocations;      // the end position and primerID of each reverse primer in the scheme
        sul::dynamic_bitset<> _ampliconOverlaps;                             // bit vector encoding all the overlap positions in the scheme
        std::unordered_map<std::string, sul::dynamic_bitset<>> _primerSites; // primer sites, stored in a bit vector per primer pool
        //std::unordered_map<std::string, Amplicon> _amplicons;           // the expected amplicons produced by the scheme
    };
} // namespace artic

#endif
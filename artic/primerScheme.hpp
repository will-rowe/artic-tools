#ifndef PRIMERSCHEME_H
#define PRIMERSCHEME_H

#include <string>
#include <unordered_map>
#include <vector>

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
        // Primer constructor
        Primer(const std::string chrom, unsigned int start, unsigned int end, const std::string primerID, const std::string poolName);

        // MergeAlt will merge a primer with an alt, yielding a primer with the maximal span
        void MergeAlt(const Primer& alt);

        // GetNumAlts returns the number of alts incorporated into this primer
        unsigned int GetNumAlts(void);

        // GetStart returns the primer start
        unsigned int GetStart(void);

        // GetEnd returns the primer end
        unsigned int GetEnd(void);

        // GetID returns the primerID
        const std::string& GetID(void) const;

        // GetDirection returns the primer direction
        signed int GetDirection(void);

        // GetBaseID returns the baseID of the primer
        const std::string& GetBaseID(void);

        // GetPrimerPool returns the pool for the primer pair
        const std::string& GetPrimerPool(void);

    private:
        std::string _chrom;
        unsigned int _start;
        unsigned int _end;
        std::string _primerID;
        std::string _poolName;
        signed int _direction; // -1 denotes reverse, 1 denotes forward
        unsigned int _numAlts; // number of alts that have been squashed into this primer object
        std::string _baseID;   // the basename of the primer (with _LEFT/_RIGHT stripped)
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
        std::pair<unsigned int, unsigned int> GetMaxSpan(void);

        // GetMinSpan returns the start and end of the amplicon, excluding the primer sequence.
        std::pair<unsigned int, unsigned int> GetMinSpan(void);

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

        // CheckScheme will check all forward primers have a paired reverse primer and record some primer scheme stats.
        void CheckScheme(void);

        // GetVersion returns the ARTIC primer scheme version (1/2/3).
        unsigned int GetVersion(void);

        // GetNumPrimers returns the total number of primers in the scheme.
        unsigned int GetNumPrimers(void);

        // GetNumAlts returns the number of alts in the primer scheme.
        unsigned int GetNumAlts(void);

        // GetNumAmplicons returns the number of amplicons in the primer scheme (after alts merged and primer pairs matched).
        unsigned int GetNumAmplicons(void);

        // GetPrimerPools returns the primer pools found in the primer scheme.
        std::vector<std::string> GetPrimerPools(void);

        // FindPrimers returns pointers to the nearest forward and reverse primer, given an alignment segment's start and end position.
        Amplicon FindPrimers(unsigned int segStart, unsigned int segEnd);

    private:
        unsigned int _version;                                               // the primer scheme version (based on the ARTIC versioning)
        unsigned int _numPrimers;                                            // the total number of primers in the scheme
        unsigned int _numAlts;                                               // the number of alts that were merged when the scheme was read
        unsigned int _numAmplicons;                                          // the number of amplicons in the scheme
        unsigned int _minStart;                                              // the minimum start value seen in the primer scheme
        unsigned int _maxEnd;                                                // the maximum end value seen in the primer scheme
        schemeMap _fPrimers;                                                 // the forward primers for the scheme
        schemeMap _rPrimers;                                                 // the reverse primers for the scheme
        std::vector<std::pair<unsigned int, std::string>> _fPrimerLocations; // the start position and primerID of each forward primer in the scheme
        std::vector<std::pair<unsigned int, std::string>> _rPrimerLocations; // the end position and primerID of each reverse primer in the scheme
        std::vector<std::string> _primerPools;                               // the primer pools found in the primer scheme
    };
}; // namespace artic

#endif
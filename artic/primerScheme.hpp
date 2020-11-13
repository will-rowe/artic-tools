#ifndef PRIMERSCHEME_H
#define PRIMERSCHEME_H

#include <boost/dynamic_bitset.hpp>
#include <htslib/faidx.h>
#include <string>
#include <vector>

#include "bytell_hash_map.hpp"
#include "kmers.hpp"

namespace artic
{
    class Primer;
    class PrimerScheme;
    class Amplicon;
    typedef ska::bytell_hash_map<std::string, Primer> primermap_t;
    typedef ska::bytell_hash_map<artic::kmer_t, std::vector<unsigned int>> kmermap_t;

    // SchemeArgs is used to pass arguments to the scheme functions.
    typedef struct SchemeArgs
    {
        std::string schemeName;
        unsigned int schemeVersion;
        std::string outDir;
        std::string refSeqFile;     // fasta file with reference sequence
        std::string schemeFile;     // bed file with primer coordinates
        std::string primerSeqsFile; // fasta file with primer sequences
        std::string insertsFile;    // bed file with insert coordinates
    } SchemeArgs;

    // DownloadScheme will download a specified primer scheme and the reference sequence.
    void DownloadScheme(SchemeArgs& args);

    // ValidateScheme will load and validate a specified primer scheme, log stats and return the primer scheme object.
    // It will optionally write primer insert coordinates and primer sequences to files.
    PrimerScheme ValidateScheme(SchemeArgs& args);

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
        Primer(unsigned int start, unsigned int end, const std::string primerID, std::size_t poolID);

        // MergeAlt will merge a primer with an alt, yielding a primer with the maximal span.
        void MergeAlt(const Primer& alt);

        // GetNumAlts returns the number of alts incorporated into this primer.
        unsigned int GetNumAlts(void) const;

        // GetStart returns the primer start.
        int64_t GetStart(void) const;

        // GetEnd returns the primer end.
        int64_t GetEnd(void) const;

        // GetLen returns the length of the primer sequence.
        unsigned int GetLen(void) const;

        // GetName returns the primerID.
        const std::string& GetName(void) const;

        // GetBaseID returns the baseID of the primer (with _LEFT/_RIGHT stripped).
        std::string GetBaseID(void) const;

        // GetPrimerPoolID returns the pool ID for the primer pair.
        std::size_t GetPrimerPoolID(void) const;

        // IsForward returns the primer direction (true = forward, false = reverse).
        bool IsForward(void);

        // GetSeq returns the primer sequence from a reference.
        void GetSeq(faidx_t* reference, const std::string& refID, std::string& primerSeq) const;

    private:
        int64_t _start;        // the reference start position of the primer (0-based, half-open)
        int64_t _end;          // the reference end position of the primer (0-based, half-open) -- NOT INCLUDED IN PRIMER
        std::string _primerID; // the full ID of the primer (_ALT will have been removed during merging)
        std::size_t _poolID;   // the primer pool ID
        bool _isForward;       // true if the forward primer, false if the reverse primer
        unsigned int _numAlts; // the number of alternate primers that have been merged into this primer
        std::size_t _baseIDit; // iterator used to shorten the primerID to the baseID (i.e. strips _LEFT/_RIGHT)
    };

    //******************************************************************************
    // PrimerScheme class handles the ARTIC style primer schemes.
    //******************************************************************************
    class PrimerScheme
    {
    public:
        // PrimerScheme constructors and destructor.
        PrimerScheme(const std::string& inputFile);
        ~PrimerScheme(void);

        // GetFileName returns the filename that the primer scheme was loaded from.
        const std::string& GetFileName(void) const;

        // GetReferenceName returns the reference sequence ID found in the primer scheme.
        const std::string& GetReferenceName(void) const;

        // GetNumPrimers returns the total number of primers in the scheme.
        unsigned int GetNumPrimers(void);

        // GetMinPrimerLen returns the minimum primer length in the scheme.
        unsigned int GetMinPrimerLen(void);

        // GetMaxPrimerLen returns the maximum primer length in the scheme.
        unsigned int GetMaxPrimerLen(void);

        // GetNumAlts returns the number of alts in the primer scheme.
        unsigned int GetNumAlts(void);

        // GetNumAmplicons returns the number of amplicons in the primer scheme (after alts merged and primer pairs matched).
        unsigned int GetNumAmplicons(void);

        // GetMeanAmpliconSpan returns the mean amplicon span (including primer sequence).
        unsigned int GetMeanAmpliconSpan(void);

        // GetMaxAmpliconSpan returns the max amplicon span (including primer sequence).
        unsigned int GetMaxAmpliconSpan(void);

        // GetPrimerPools returns all the primer pools found in the primer scheme.
        std::vector<std::string> GetPrimerPools(void);

        // GetPrimerPool returns the primer pool for the provided pool ID.
        const std::string& GetPrimerPool(std::size_t poolID) const;

        // GetPrimerPoolID returns the primer pool ID for the provided pool name.
        std::size_t GetPrimerPoolID(const std::string& poolName) const;

        // GetRefStart returns the first position in the reference covered by the primer scheme.
        int64_t GetRefStart(void);

        // GetRefEnd returns the last position in the reference covered by the primer scheme.
        int64_t GetRefEnd(void);

        // GetNumOverlaps returns the number of reference positions covered by more than one amplicon.
        unsigned int GetNumOverlaps(void);

        // GetExpAmplicons returns a vector containing the amplicons the scheme expects to produce.
        const std::vector<Amplicon>& GetExpAmplicons(void);

        // GetAmpliconName returns a string name for an amplicon in the scheme, based on the provided amplicon int ID.
        const std::string GetAmpliconName(unsigned int id);

        // FindPrimers returns pointers to the nearest forward and reverse primer, given an alignment segment's start and end position.
        Amplicon FindPrimers(int64_t segStart, int64_t segEnd);

        // CheckAmpliconOverlap returns true if the queried position is covered by multiple amplicons (incl. primer sequence).
        bool CheckAmpliconOverlap(int64_t pos);

        // CheckPrimerSite returns true if the queried position is a primer site for the given pool.
        bool CheckPrimerSite(int64_t pos, const std::string& poolName);

        // GetPrimerKmers will int encode k-mers from all primers in the scheme and deposit them in the provided map, linked to their amplicon primer origin(s).
        void GetPrimerKmers(const std::string& reference, uint32_t kSize, kmermap_t& kmerMap);

    private:
        void _loadScheme(const std::string& filename);                  // _loadScheme will load an input file and create the primer objects.
        void _validateScheme(void);                                     // _validateScheme will check all forward primers have a paired reverse primer and record some primer scheme stats.
        std::string _filename;                                          // the file that the scheme was loaded from
        std::string _referenceID;                                       // the ID of the reference sequence covered by the primer scheme
        unsigned int _numPrimers;                                       // the total number of primers in the scheme
        unsigned int _numAlts;                                          // the number of alts that were merged when the scheme was read
        unsigned int _numAmplicons;                                     // the number of amplicons in the scheme
        unsigned int _meanAmpliconSpan;                                 // the mean amplicon span (incl. primers)
        unsigned int _maxAmpliconSpan;                                  // the max amplicon span (incl. primers)
        unsigned int _minPrimerLen;                                     // the minimum primer length in the scheme
        unsigned int _maxPrimerLen;                                     // the maximum primer length in the scheme
        int64_t _refStart;                                              // the first position in the reference covered by the primer scheme
        int64_t _refEnd;                                                // the last position in the reference covered by the primer scheme
        std::vector<std::string> _primerPools;                          // the primer pool IDs found in the primer scheme
        primermap_t _fPrimers;                                          // the forward primers for the scheme
        primermap_t _rPrimers;                                          // the reverse primers for the scheme
        std::vector<std::pair<int64_t, std::string>> _fPrimerLocations; // the start position and primerID of each forward primer in the scheme
        std::vector<std::pair<int64_t, std::string>> _rPrimerLocations; // the end position and primerID of each reverse primer in the scheme
        boost::dynamic_bitset<> _ampliconOverlaps;                      // bit vector encoding all the overlap positions in the scheme
        boost::dynamic_bitset<> _primerSites;                           // primer sites, stored in a bit vector and offset by primer pool ID
        std::vector<Amplicon> _expAmplicons;                            // the expected amplicons produced by the scheme
    };

    //******************************************************************************
    // Amplicon is a container for two primers.
    //******************************************************************************
    class Amplicon
    {
    public:
        // Amplicon constructor.
        Amplicon(Primer* p1, Primer* p2);

        // SetID will assign an ID to the amplicon.
        void SetID(unsigned int id);

        // IsProperlyPaired returns true if the amplicon primers are properly paired.
        bool IsProperlyPaired(void);

        // GetName returns the name for the amplicon (combines primer IDs).
        const std::string GetName(void) const;

        // GetID returns the numberical ID for the amplicon.
        unsigned int GetID(void) const;

        // GetPrimerPool returns the pool ID for the primer pair (0 returned if not properly paired).
        std::size_t GetPrimerPoolID(void);

        // GetMaxSpan returns the start and end of the amplicon, including the primer sequence.
        std::pair<int64_t, int64_t> GetMaxSpan(void);

        // GetMinSpan returns the start and end of the amplicon, excluding the primer sequence.
        std::pair<int64_t, int64_t> GetMinSpan(void);

        // GetForwardPrimer returns a pointer to the forward primer in the amplicon.
        const Primer* GetForwardPrimer(void);

        // GetReversePrimer returns a pointer to the reverse primer in the amplicon.
        const Primer* GetReversePrimer(void);

        // AddKmers adds the k-mers from a sequence to the amplicon.
        void AddKmers(const char* seq, uint32_t seqLen, uint32_t kSize);

    private:
        Primer* _fPrimer;        // pointer to the forward primer object
        Primer* _rPrimer;        // pointer to the reverse primer object
        bool _isProperlyPaired;  // denotes if amplicon has properly paired primers
        unsigned int _id;        // the amplicon identifier
        artic::kmerset_t _kmers; // the set of k-mers for this amplicon
    };

} // namespace artic

#endif
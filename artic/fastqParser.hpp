#ifndef FASTQ_PARSER_H
#define FASTQ_PARSER_H

#include <vector>
#include <zlib.h>

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);
#endif

namespace artic
{
    class FastqReader
    {

    public:
        FastqReader(const std::vector<std::string> fnames);
        ~FastqReader();
        void Close();
        void Reopen();
        int GetRecord(std::string& seq, size_t& id);

    private:
        std::vector<std::string>::const_iterator openNextFile();
        std::vector<std::string> _files;                   // the FASTQ filenames to read
        std::vector<std::string>::const_iterator _fileItr; // an iterator for the filenames
        unsigned int _fileNumber;                          // the number of files read
        gzFile _file;                                      // the current open file
        kseq_t* _kseq;
    };
} // namespace artic

#endif // FASTQ_PARSER_H
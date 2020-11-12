#include <iostream>
#include <sstream>
#include <string>

#include "fastqParser.hpp"

//
artic::FastqReader::FastqReader(const std::vector<std::string> files)
    : _files(files), _fileNumber(0), _kseq(NULL)
{
    _fileItr = files.begin();
    _file = gzopen(_fileItr->c_str(), "r");
    _kseq = kseq_init(_file);
}

//
artic::FastqReader::~FastqReader()
{
    closeFile();
    _fileItr = _files.end();
}

// returns >=0 (length of seq), -1 end of last file, -2 truncated quality string
int artic::FastqReader::GetRecord(std::string& seq, size_t& id)
{
    const int r = kseq_read(_kseq);
    if (r >= 0)
        seq.assign(_kseq->seq.s);
    else if (r == -1)
    {
        if (openNextFile())
        {
            id = _fileNumber;
            return GetRecord(seq, id);
        }
    }
    return r;
}

//
bool artic::FastqReader::openNextFile()
{
    // close current file
    closeFile();

    // return false if there are no more files to open
    if (_fileNumber == _files.size())
        return false;

    // otherwise open another file
    _file = gzopen(_fileItr->c_str(), "r");
    _kseq = kseq_init(_file);
    return true;
}

//
void artic::FastqReader::closeFile()
{
    if (_kseq == NULL)
        return;
    kseq_destroy(_kseq);
    gzclose(_file);
    _fileItr++;
    _fileNumber++;
    _kseq = NULL;
    return;
}
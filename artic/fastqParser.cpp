#include <iostream>
#include <sstream>
#include <string>

#include "fastqParser.hpp"

//
FastqFile::FastqFile(const std::vector<std::string> files)
    : _files(files), _fileNumber(0), _kseq(NULL)
{
    _fileItr = files.begin();
    _file = gzopen(_fileItr->c_str(), "r");
    _kseq = kseq_init(_file);
}

//
FastqFile::~FastqFile()
{
    Close();
}

//
void FastqFile::Close()
{
    if (_kseq != NULL)
    {

        kseq_destroy(_kseq);
        gzclose(_file);
        _fileItr = _files.end();
        _kseq = NULL;
    }
}

//
void FastqFile::Reopen()
{
    Close();
    _fileItr = _files.begin();
    _file = gzopen(_fileItr->c_str(), "r");
    _kseq = kseq_init(_file);
}

// returns >=0 (length of seq), -1 end of last file, -2 truncated quality string
int FastqFile::GetRecord(std::string& seq, size_t& id)
{
    const int r = kseq_read(_kseq);

    if (r >= 0)
        seq.assign(_kseq->seq.s);
    else if (r == -1)
    {

        openNextFile();

        if (_fileItr != _files.end())
        {

            id = _fileNumber;

            return GetRecord(seq, id);
        }
    }

    return r;
}

//
std::vector<std::string>::const_iterator FastqFile::openNextFile()
{
    if (_fileItr != _files.end())
    {
        kseq_destroy(_kseq);
        gzclose(_file);
        _kseq = NULL;

        ++_fileItr;
        ++_fileNumber;
        if (_fileItr != _files.end())
        {
            _file = gzopen(_fileItr->c_str(), "r");
            _kseq = kseq_init(_file);
        }
    }
    return _fileItr;
}
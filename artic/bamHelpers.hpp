#ifndef BAMUTILS_H
#define BAMUTILS_H

#include <htslib/sam.h>
#include <string>
#include <vector>

namespace artic
{
    // TrimAlignment will softmask an alignment from its start/end up to the provided position.
    void TrimAlignment(bam1_t* record, unsigned int maskEnd, bool reverse);

    // AddPGtoHeader will add a PG tag to an existing BAM header based on the provided command.
    void AddPGtoHeader(bam_hdr_t** header, const std::string& userCmd);

    // AddRGtoHeader will add a RG tag to an existing BAM header based on the provided read group.
    void AddRGtoHeader(bam_hdr_t** header, std::string& rg);

} // namespace artic

#endif
#ifndef VERSION_H
#define VERSION_H

#define VERSION_MAJOR 0
#define VERSION_MINOR 1
#define VERSION_PATCH 0

#include <string>

const std::string PROG_NAME = "artic_tools";

namespace artic
{
    std::string GetVersion();
}; // namespace artic

#endif

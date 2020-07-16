#include "version.hpp"

// GetVersion will return a formatted version string.
std::string artic::GetVersion()
{
    char buffer[5];
    std::sprintf(buffer, "%d.%d.%d", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH);
    return buffer;
};

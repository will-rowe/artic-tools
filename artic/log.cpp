#include "spdlog/sinks/stdout_color_sinks.h"

#include "log.hpp"
#include "version.hpp"

std::shared_ptr<spdlog::logger> artic::Log::s_Corelogger;
std::shared_ptr<spdlog::logger> artic::Log::s_Clientlogger;

void artic::Log::Init(std::string subtool)
{
    spdlog::set_pattern("[%H:%M:%S] [%n] %v");

    s_Corelogger = spdlog::stderr_color_mt("CORE");
    s_Corelogger->set_level(spdlog::level::trace);
    std::string logcaller = PROG_NAME + "::" + subtool;
    s_Clientlogger = spdlog::stderr_color_mt(logcaller);
    s_Clientlogger->set_level(spdlog::level::trace);
}

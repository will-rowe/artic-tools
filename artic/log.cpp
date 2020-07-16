#include "spdlog/sinks/stdout_color_sinks.h"

#include "log.hpp"
#include "version.hpp"

std::shared_ptr<spdlog::logger> artic::Log::s_Corelogger;
std::shared_ptr<spdlog::logger> artic::Log::s_Clientlogger;

void artic::Log::Init()
{
    spdlog::set_pattern("%^[%T] %n: %v%$");

    s_Corelogger = spdlog::stderr_color_mt("CORE");
    s_Corelogger->set_level(spdlog::level::trace);

    s_Clientlogger = spdlog::stderr_color_mt(PROG_NAME);
    s_Clientlogger->set_level(spdlog::level::trace);
}

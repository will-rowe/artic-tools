#ifndef LOGGER_H
#define LOGGER_H

#include "spdlog/spdlog.h"
#include <memory>

namespace artic
{
    class Log
    {
    private:
        static std::shared_ptr<spdlog::logger> s_Corelogger;
        static std::shared_ptr<spdlog::logger> s_Clientlogger;

    public:
        static void Init(std::string subtool);
        inline static std::shared_ptr<spdlog::logger>& GetCoreLogger() { return s_Corelogger; };
        inline static std::shared_ptr<spdlog::logger>& GetClientLogger() { return s_Clientlogger; };
    };

// Core Logger Macro
#define LOG_CORE_TRACE(...) artic::Log::GetCoreLogger()->trace(__VA_ARGS__);
#define LOG_CORE_INFO(...) artic::Log::GetCoreLogger()->info(__VA_ARGS__);
#define LOG_CORE_WARN(...) artic::Log::GetCoreLogger()->warn(__VA_ARGS__);
#define LOG_CORE_ERROR(...) artic::Log::GetCoreLogger()->error(__VA_ARGS__);
#define LOG_CORE_FATAL(...) artic::Log::GetCoreLogger()->fatal(__VA_ARGS__);

// Client Logger Macro
#define LOG_TRACE(...) artic::Log::GetClientLogger()->trace(__VA_ARGS__);
#define LOG_INFO(...) artic::Log::GetClientLogger()->info(__VA_ARGS__);
#define LOG_WARN(...) artic::Log::GetClientLogger()->warn(__VA_ARGS__);
#define LOG_ERROR(...) artic::Log::GetClientLogger()->error(__VA_ARGS__);
#define LOG_FATAL(...) artic::Log::GetClientLogger()->fatal(__VA_ARGS__);

} // namespace artic

#endif
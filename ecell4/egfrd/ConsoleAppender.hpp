#ifndef ECELL4_EGFRD_CONSOLE_LOGGER_HPP
#define ECELL4_EGFRD_CONSOLE_LOGGER_HPP
#include "Logger.hpp"
#include <cstdio>
#include <cstring>

namespace ecell4
{
namespace egfrd
{

class ConsoleAppender: public LogAppender
{
public:
    typedef LogAppender base_type;

public:
    ~ConsoleAppender() override = default;

    void flush() override
    {
        std::fflush(stderr);
    }

    void operator()(enum Logger::level lv, char const* name, char const** chunks)
    {
        std::fprintf(stderr, "%s: %-8s ", name, Logger::stringize_error_level(lv));
        for (char const** p = chunks; *p; ++p)
        {
            std::fwrite(*p, sizeof(char), std::strlen(*p), stderr);
        }
        std::fputc('\n', stderr);
    }
};

} // egfrd
} // ecell4
#endif /* CONSOLE_LOGGER_HPP */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <cstdio>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "ConsoleLogger.hpp"

ConsoleLogger::~ConsoleLogger() {}

void ConsoleLogger::logv(enum level lv, char const* format, va_list ap)
{
    using namespace boost::posix_time;
    std::fprintf(stderr, "[%s] %s: %-8s ", to_iso_string(second_clock::local_time()).c_str(), name_.c_str(), stringize_error_level(lv));
    std::vfprintf(stderr, format, ap);
    std::fputc('\n', stderr);
}

char const* ConsoleLogger::stringize_error_level(enum level lv)
{
    static char const* names[] = {
        "OFF",
        "DEBUG",
        "INFO",
        "WARN",
        "ERROR",
        "FATAL"
    };
    return static_cast<std::size_t>(lv) >= sizeof(names) / sizeof(*names) ? "???": names[lv];
}

ConsoleLoggerFactory::~ConsoleLoggerFactory() {}

Logger* ConsoleLoggerFactory::operator()(char const* logger_name) const
{
    return new ConsoleLogger(logger_name);
}

char const* ConsoleLoggerFactory::get_name() const
{
    return "ConsoleLogger";
}

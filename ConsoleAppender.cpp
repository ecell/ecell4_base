#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <cstdio>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "ConsoleAppender.hpp"

ConsoleAppender::~ConsoleAppender() {}

void ConsoleAppender::operator()(enum Logger::level lv, boost::posix_time::ptime const& tm, char const* name, char const** chunks)
{
    using namespace boost::posix_time;
    std::fprintf(stderr, "[%s] %s: %-8s ", to_iso_string(tm).c_str(),
                 name, Logger::stringize_error_level(lv));
    for (char const** p = chunks; *p; ++p)
        std::fwrite(*p, sizeof(char), strlen(*p), stderr);
    std::fputc('\n', stderr);
}

void ConsoleAppender::flush()
{
    std::fflush(stderr);
}

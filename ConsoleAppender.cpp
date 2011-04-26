#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <cstdio>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "ConsoleAppender.hpp"

ConsoleAppender::~ConsoleAppender() {}

void ConsoleAppender::operator()(enum Logger::level lv, char const* name, char const** chunks)
{
    std::fprintf(stderr, "%s: %-8s ",
      name, Logger::stringize_error_level(lv));
    for (char const** p = chunks; *p; ++p)
        std::fwrite(*p, sizeof(char), strlen(*p), stderr);
    std::fputc('\n', stderr);
}

void ConsoleAppender::flush()
{
    std::fflush(stderr);
}

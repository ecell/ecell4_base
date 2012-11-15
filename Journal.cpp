#include <cstdio>
#include <cstring>

#include "Journal.hpp"


namespace ecell4
{

Journal::Journal(char const* name)
    : name_(name)
{
    ;
}

Journal::~Journal()
{
    ;
}

void Journal::level(enum Journal::level level)
{
    ensure_initialized();
    level_ = level;
}

enum Journal::level Journal::level() const
{
    const_cast<Journal*>(this)->ensure_initialized();
    return level_;
}

void Journal::logv(enum level lv, char const* format, va_list ap)
{
    ensure_initialized();

    if (lv < level_)
        return;

    char buf[1024];
    std::vsnprintf(buf, sizeof(buf), format, ap);

    std::fprintf(stderr, "%s: %-8s ", name_.c_str(), stringize_error_level(lv));
    std::fwrite(buf, sizeof(char), std::strlen(buf), stderr);
    std::fputc('\n', stderr);
}

void Journal::flush()
{
    ensure_initialized();

    std::fflush(stderr);
}

char const* Journal::stringize_error_level(enum level lv)
{
    static char const* names[] = {
        "OFF",
        "DEBUG",
        "INFO",
        "WARN",
        "ERROR",
        "FATAL"
    };

    return (static_cast<std::size_t>(lv) >= sizeof(names) / sizeof(*names)
            ? "???": names[lv]);
}


// Journal& Journal::get_journal(char const* name)
// {
//     ;
// }

inline void Journal::ensure_initialized()
{
    ; // not implemented yet
}

} // ecell4

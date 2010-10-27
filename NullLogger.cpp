#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "NullLogger.hpp"

NullLogger::~NullLogger() {}

void NullLogger::logv(enum level lv, char const* format, va_list ap)
{
}

void NullLogger::flush()
{
}

NullLoggerFactory::~NullLoggerFactory() {}

Logger* NullLoggerFactory::operator()(char const* logger_name) const
{
    return new NullLogger(logger_name);
}

char const* NullLoggerFactory::get_name() const
{
    return "NullLogger";
}

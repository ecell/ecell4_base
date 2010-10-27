#ifndef NULL_LOGGER_HPP
#define NULL_LOGGER_HPP

#include <string>

#include "Logger.hpp"

class NullLogger: public Logger
{
public:
    virtual ~NullLogger();

    virtual void logv(enum level lv, char const* format, va_list ap);

    virtual void flush();

    NullLogger(char const* name): name_(name) {}

private:
    static char const* stringize_error_level(enum level lv);

private:
    std::string name_;
};

class NullLoggerFactory: public LoggerFactory
{
public:
    virtual ~NullLoggerFactory();

    virtual Logger* operator()(char const* logger_name) const;

    virtual char const* get_name() const;
};

#endif /* NULL_LOGGER_HPP */

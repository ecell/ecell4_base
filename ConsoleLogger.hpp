#ifndef CONSOLE_LOGGER_HPP
#define CONSOLE_LOGGER_HPP

#include <string>

#include "Logger.hpp"

class ConsoleLogger: public Logger
{
public:
    virtual ~ConsoleLogger();

    virtual void logv(enum level lv, char const* format, va_list ap);

    ConsoleLogger(char const* name): name_(name) {}

private:
    static char const* stringize_error_level(enum level lv);

private:
    std::string name_;
};

class ConsoleLoggerFactory: public LoggerFactory
{
public:
    virtual ~ConsoleLoggerFactory();

    virtual Logger* operator()(char const* logger_name) const;

    virtual char const* get_name() const;
};

#endif /* CONSOLE_LOGGER_HPP */

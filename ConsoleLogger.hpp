#ifndef CONSOLE_LOGGER_HPP
#define CONSOLE_LOGGER_HPP

#include <string>

#include "Logger.hpp"

class ConsoleLogger: public Logger
{
public:
    virtual ~ConsoleLogger();

    virtual void level(enum level);

    virtual enum level level() const;

    virtual void logv(enum level lv, char const* format, va_list ap);

    virtual void flush();

    ConsoleLogger(char const* name);

private:
    static char const* stringize_error_level(enum level lv);

private:
    std::string name_;
    enum level level_;
};

class ConsoleLoggerFactory: public LoggerFactory
{
public:
    virtual ~ConsoleLoggerFactory();

    virtual void level(enum Logger::level level);

    virtual Logger* operator()(char const* logger_name) const;

    virtual char const* get_name() const;

    ConsoleLoggerFactory();

private:
    enum Logger::level level_;
};

#endif /* CONSOLE_LOGGER_HPP */

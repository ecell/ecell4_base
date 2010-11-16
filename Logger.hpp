#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <cstdarg>
#include <boost/shared_ptr.hpp>

class Logger
{
public:
    enum level
    {
        L_OFF = 0,
        L_DEBUG = 1,
        L_INFO = 2,
        L_WARNING = 3,
        L_ERROR = 4,
        L_FATAL = 5
    };

public:
    virtual ~Logger();

    void debug(char const* format, ...)
    {
        va_list ap;
        va_start(ap, format);
        logv(L_DEBUG, format, ap);
        va_end(ap);
    }

    void info(char const* format, ...)
    {
        va_list ap;
        va_start(ap, format);
        logv(L_INFO, format, ap);
        va_end(ap);
    }

    void warn(char const* format, ...)
    {
        va_list ap;
        va_start(ap, format);
        logv(L_WARNING, format, ap);
        va_end(ap);
    }

    void error(char const* format, ...)
    {
        va_list ap;
        va_start(ap, format);
        logv(L_ERROR, format, ap);
        va_end(ap);
    }

    void fatal(char const* format, ...)
    {
        va_list ap;
        va_start(ap, format);
        logv(L_FATAL, format, ap);
        va_end(ap);
    }

    void log(enum level lv, char const* format, ...)
    {
        va_list ap;
        va_start(ap, format);
        logv(lv, format, ap);
        va_end(ap);
    }

    virtual void level(enum level) = 0;

    virtual enum level level() const = 0;

    virtual void logv(enum level lv, char const* format, va_list ap) = 0;

    virtual void flush() = 0;

    static Logger& get_logger(char const* name);
};

class LoggerFactory
{
public:
    virtual ~LoggerFactory();

    virtual void level(enum Logger::level level) = 0;

    virtual Logger* operator()(char const* logger_name) const = 0;

    virtual char const* get_name() const = 0;

    static void register_logger_factory(
            char const* logger_name_pattern,
            boost::shared_ptr<LoggerFactory> const&);

    static boost::shared_ptr<LoggerFactory>
           get_logger_factory(char const* logger_name);
};

#ifdef DEBUG
#   define LOG_DEBUG(args) if (log_.level() == Logger::L_DEBUG) log_.debug args
#else
#   define LOG_DEBUG(args)
#endif 

#define LOG_INFO(args) if (enum Logger::level const level = log_.level()) if (level <= Logger::L_INFO) log_.info args

#endif /* LOGGER_HPP */

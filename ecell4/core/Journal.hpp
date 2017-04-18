#ifndef ECELL4_JOURNAL_HPP
#define ECELL4_JOURNAL_HPP

#include <cstdarg>
#include <string>


namespace ecell4
{

class Journal
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

    Journal(char const* name);

    ~Journal();

    void level(enum level level);
    enum level level() const;

    char const* name() const
    {
        return name_.c_str();
    }

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

    void logv(enum level lv, char const* format, va_list ap);
    void flush();

    static char const* stringize_error_level(enum level lv);

    // static Journal& get_journal(char const* name);

private:

    void ensure_initialized();

protected:

    const std::string name_;
    enum level level_;
};

} // ecell4

#endif /* ECELL4_JOURNAL_HPP */

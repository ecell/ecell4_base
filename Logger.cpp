#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <string>
#include <map>
#include <cstdio>
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "Logger.hpp"
#include "utils/pair.hpp"
#include "utils/fun_composition.hpp"
#include "utils/fun_wrappers.hpp"
#include "assoc_container_traits.hpp"
#include "map_adapter.hpp"

struct map_adapter_handler
{
    template<typename Tadapter_>
    void destroy(Tadapter_& cntnr) const
    {
        std::for_each(boost::begin(cntnr), boost::end(cntnr),
                compose_unary(
                    delete_ptr<
                        typename boost::remove_pointer<
                            typename Tadapter_::mapped_type>::type>(),
                    select_second<typename Tadapter_::value_type>()));
    }

    template<typename Tadapter_, typename Titer_>
    void insert(Titer_ const& b, Titer_ const& e) const
    {
    }

    template<typename Tadapter_>
    void insert(typename Tadapter_::value_type const& val) const
    {
    }
};

Logger::~Logger()
{
}

Logger& Logger::get_logger(char const* name)
{
    typedef map_adapter<std::map<std::string, Logger*>, map_adapter_handler> loggers_type;
    static map_adapter_handler hdlr;
    static loggers_type loggers(hdlr);
    std::string _name(name);
    std::pair<loggers_type::iterator, bool> i(
            loggers.insert(loggers_type::value_type(_name, 0)));

    if (i.second)
    {
        Logger* log = LoggerFactory::get_logger_factory(name).create();
        log->set_name(name);
        (*i.first).second = log;
    }

    return *(*i.first).second;
}

LoggerFactory::~LoggerFactory()
{
}

class ConsoleLogger: public Logger
{
public:
    virtual ~ConsoleLogger() {}

    virtual void set_name(char const* name)
    {
        name_ = name;
    }

    virtual void logv(enum level lv, char const* format, va_list ap)
    {
        using namespace boost::posix_time;
        std::fprintf(stderr, "[%s] %s: %-8s ", to_iso_string(second_clock::local_time()).c_str(), name_.c_str(), stringize_error_level(lv));
        std::vfprintf(stderr, format, ap);
        std::fputc('\n', stderr);
    }

private:
    static char const* stringize_error_level(enum level lv)
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

private:
    std::string name_;
};

class ConsoleLoggerFactory: public LoggerFactory
{
public:
    virtual ~ConsoleLoggerFactory() {}

    virtual Logger* create() const
    {
        return new ConsoleLogger();
    }
};

LoggerFactory& LoggerFactory::get_logger_factory(char const* name)
{
    static ConsoleLoggerFactory singleton_;
    return singleton_;
}

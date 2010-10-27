#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <string>
#include <map>
#include <utility>
#include <boost/regex.hpp>
#include <boost/foreach.hpp>
#include <boost/type_traits/remove_pointer.hpp>
#include "Logger.hpp"
#include "ConsoleLogger.hpp"
#include "utils/pair.hpp"
#include "utils/fun_composition.hpp"
#include "utils/fun_wrappers.hpp"
#include "utils/assoc_container_traits.hpp"
#include "utils/map_adapter.hpp"

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

struct LoggerProxy: public Logger
{
    virtual ~LoggerProxy() {}

    virtual void logv(enum level lv, char const* format, va_list ap)
    {
        fetch_logger();
        logger_->logv(lv, format, ap);
    }

    virtual void flush()
    {
        fetch_logger();
        logger_->flush();
    }

    void fetch_logger()
    {
        if (logger_)
            return;
        logger_ = (*factory_)(name_);
    }

    LoggerProxy(boost::shared_ptr<LoggerFactory> factory, char const* name)
        : factory_(factory), logger_(0), name_(name) {}

private:
    boost::shared_ptr<LoggerFactory> factory_;
    Logger* logger_;
    char const* name_;
};

struct DeferredLoggerFactory: public LoggerFactory
{
    virtual ~DeferredLoggerFactory() {}

    virtual Logger* operator()(char const* name) const
    {
        const_cast<DeferredLoggerFactory*>(this)->fetch_factory();
        Logger* const log((*factory_)(name));
        return log;
    }

    virtual char const* get_name() const
    {
        const_cast<DeferredLoggerFactory*>(this)->fetch_factory();
        return factory_->get_name();
    }

    DeferredLoggerFactory(char const* name): name_(name), factory_() {}

private:
    void fetch_factory()
    {
        if (!factory_)
            factory_ = get_logger_factory(name_);
        BOOST_ASSERT(factory_.get());
    }

private:
    char const* name_;
    boost::shared_ptr<LoggerFactory> factory_;
};

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
        Logger* log = new LoggerProxy(
            boost::shared_ptr<LoggerFactory>(
                new DeferredLoggerFactory(name)), name);
        (*i.first).second = log;
    }

    return *(*i.first).second;
}

LoggerFactory::~LoggerFactory()
{
}

class LoggerFactoryRegistry
{
private:
    typedef std::pair<boost::regex, boost::shared_ptr<LoggerFactory> > entry_type;
public:
    void register_logger_factory(char const* logger_name_pattern,
                                 boost::shared_ptr<LoggerFactory> const& factory)
    {
        factories_.push_back(entry_type(boost::regex(logger_name_pattern), factory));
    }

    boost::shared_ptr<LoggerFactory>
    get_logger_factory(char const* logger_name) const
    {
        BOOST_FOREACH (entry_type const& i, factories_)
        {
            if (boost::regex_match(logger_name,
                                   logger_name + std::strlen(logger_name),
                                   i.first))
            {
                return i.second;
            }
        }
        BOOST_ASSERT(default_factory_.get());
        return default_factory_;
    }

    LoggerFactoryRegistry(): default_factory_(new ConsoleLoggerFactory()) {}

private:
    std::vector<entry_type> factories_;
    boost::shared_ptr<ConsoleLoggerFactory> default_factory_;
};

static LoggerFactoryRegistry registry;
    
void LoggerFactory::register_logger_factory(
        char const* logger_name_pattern,
        boost::shared_ptr<LoggerFactory> const& factory)
{
    registry.register_logger_factory(logger_name_pattern, factory);
}


boost::shared_ptr<LoggerFactory>
LoggerFactory::get_logger_factory(char const* name)
{
    return registry.get_logger_factory(name);
}

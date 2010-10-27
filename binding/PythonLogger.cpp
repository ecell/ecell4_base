#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <map>
#include <cstdio>
#include "PythonLogger.hpp"
#include "binding_common.hpp"

namespace binding {

static boost::python::object logging_module;
std::map<enum Logger::level, boost::python::object> loglevel_map;

class PythonLogger: public Logger
{
public:
    virtual ~PythonLogger() {}

    virtual void logv(enum level lv, char const* format, va_list ap)
    {
        char buf[2048];
        std::vsnprintf(buf, sizeof(buf), format, ap);
        logger_.attr("log")(loglevel_map[lv], "%s", buf);
    }

    PythonLogger(boost::python::object logger): logger_(logger) {}

private:
    boost::python::object logger_;
};

class PythonLoggerFactory: public LoggerFactory
{
public:
    virtual ~PythonLoggerFactory();

    virtual Logger* operator()(char const* name) const;

    virtual char const* get_name() const;
};

PythonLoggerFactory::~PythonLoggerFactory()
{
}
Logger* PythonLoggerFactory::operator()(char const* name) const
{
    using namespace boost::python;
    object getLogger(getattr(logging_module, "getLogger"));
    return new PythonLogger(getLogger(name));
}

char const* PythonLoggerFactory::get_name() const
{
    return "PythonLogger";
}

boost::python::objects::class_base
register_logger_factory_class(char const* name)
{
    using namespace boost::python;
    typedef LoggerFactory impl_type;

    return class_<impl_type, boost::shared_ptr<impl_type>, boost::noncopyable>(name, no_init)
        .add_property("name", &impl_type::get_name)
        .def("register_logger_factory", &impl_type::register_logger_factory)
        .staticmethod("register_logger_factory")
        ;
}

boost::python::objects::class_base
register_python_logger_factory_class(char const* name)
{
    using namespace boost::python;
    typedef PythonLoggerFactory impl_type;

    logging_module = boost::python::object(boost::python::borrowed(PyImport_Import(PyString_FromString("logging"))));
    if (PyErr_Occurred())
    {
        throw_error_already_set();
    }

    loglevel_map[Logger::L_OFF]     = getattr(logging_module, "NOTSET");
    loglevel_map[Logger::L_DEBUG]   = getattr(logging_module, "DEBUG");
    loglevel_map[Logger::L_INFO]    = getattr(logging_module, "INFO");
    loglevel_map[Logger::L_WARNING] = getattr(logging_module, "WARNING");
    loglevel_map[Logger::L_ERROR]   = getattr(logging_module, "ERROR");
    loglevel_map[Logger::L_FATAL]   = getattr(logging_module, "CRITICAL");

    return class_<impl_type, bases<LoggerFactory>, boost::shared_ptr<impl_type>, boost::noncopyable>(name)
        ;
}

} // namespace binding

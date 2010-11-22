#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/python.hpp>
#include "peer/converters/sequence.hpp"
#include "binding_common.hpp"

namespace binding {

boost::python::objects::enum_base
register_logger_level_enum(char const* name)
{
    using namespace boost::python;
    return enum_<enum   Logger::level>(name)
        .value("OFF",     Logger::L_OFF)
        .value("DEBUG",   Logger::L_DEBUG)
        .value("INFO",    Logger::L_INFO)
        .value("WARNING", Logger::L_WARNING)
        .value("ERROR",   Logger::L_ERROR)
        .value("FATAL",   Logger::L_FATAL)
        ;
}

boost::python::objects::class_base
register_logger_manager_class(char const* name)
{
    using namespace boost::python;
    typedef LoggerManager impl_type;

    peer::converters::register_range_to_tuple_converter<
        std::vector<boost::shared_ptr<LogAppender> > >();

    return class_<impl_type, boost::shared_ptr<impl_type>, boost::noncopyable>(name, no_init)
        .add_property("level",
            static_cast<enum Logger::level(impl_type::*)() const>(&impl_type::level),
            static_cast<void(impl_type::*)(enum Logger::level)>(&impl_type::level))
        .add_property("name", &impl_type::name)
        .add_property("appenders",
            make_function(&impl_type::appenders,
                return_value_policy<return_by_value>()))
        .def("add_appender", &impl_type::add_appender)
        .def("register_logger_manager", &impl_type::register_logger_manager)
        .staticmethod("register_logger_manager")
        .def("get_logger_manager", &impl_type::get_logger_manager)
        .staticmethod("get_logger_manager")
        ;
}

} // namespace binding

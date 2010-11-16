#ifndef PYTHON_LOGGER_HPP
#define PYTHON_LOGGER_HPP

#include <boost/python.hpp>

namespace binding {

boost::python::objects::class_base
register_logger_factory_class(char const* name);

boost::python::object
register_logger_handler_class(char const* name);

boost::python::objects::class_base
register_python_logger_factory_class(char const* name);

boost::python::objects::class_base
register_null_logger_factory_class(char const* name);

boost::python::objects::enum_base
register_logger_level_enum(char const* name);

} // namespace binding

#endif /* PYTHON_LOGGER_HPP */

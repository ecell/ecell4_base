#ifndef PYTHON_APPENDER_HPP
#define PYTHON_APPENDER_HPP

#include <boost/python.hpp>

namespace binding {

boost::python::object
register_logger_handler_class(char const* name);

boost::python::objects::class_base
register_python_appender_class(char const* name);

} // namespace binding

#endif /* PYTHON_APPENDER_HPP */

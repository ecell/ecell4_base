#ifndef CONSOLE_APPENDER_HPP
#define CONSOLE_APPENDER_HPP

#include <boost/python.hpp>

namespace binding {

boost::python::objects::class_base
register_console_appender_class(char const* name);

} // namespace binding

#endif /* CONSOLE_APPENDER_HPP */

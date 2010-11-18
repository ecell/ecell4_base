#ifndef LOG_APPENDER_HPP
#define LOG_APPENDER_HPP

#include <boost/python.hpp>

namespace binding {

boost::python::objects::class_base
register_log_appender_class(char const* name);

} // namespace binding

#endif /* LOG_APPENDER_HPP */

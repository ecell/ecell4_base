#ifndef LOGGER_MANAGER_HPP
#define LOGGER_MANAGER_HPP

#include <boost/python.hpp>

namespace binding {

boost::python::objects::class_base
register_logger_manager_class(char const* name);

boost::python::objects::enum_base
register_logger_level_enum(char const* name);

} // namespace binding

#endif /* LOGGER_MANAGER_HPP */

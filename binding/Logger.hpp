#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <boost/python.hpp>

namespace binding {

boost::python::objects::class_base
register_logger_class(char const* name);

} // namespace binding

#endif /* LOGGER_HPP */

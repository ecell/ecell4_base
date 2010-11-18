#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/python.hpp>
#include "peer/util/shared_const_ptr.hpp"
#include "binding_common.hpp"

namespace binding {

boost::python::objects::class_base
register_logger_class(char const* name)
{
    using namespace boost::python;
    typedef Logger impl_type;

    peer::util::register_shared_const_ptr_from_python<LoggerManager>();
    peer::util::register_shared_const_ptr_to_python<LoggerManager>();

    return class_<impl_type, boost::noncopyable>(name, no_init)
        .add_property("level",
            static_cast<enum impl_type::level(impl_type::*)() const>(
                &impl_type::level),
            static_cast<void (impl_type::*)(enum impl_type::level)>(
                &impl_type::level))
        .add_property("name", &impl_type::name)
        .add_property("manager",
            make_function(&impl_type::manager,
                return_value_policy<return_by_value>()))
        .def("flush", &impl_type::flush)
        .def("get_logger", &impl_type::get_logger,
            return_value_policy<reference_existing_object>())
        .staticmethod("get_logger")
        .def("stringize_error_level", &impl_type::stringize_error_level)
        .staticmethod("stringize_error_level")
        ;
}

} // namespace binding

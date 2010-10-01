#ifndef BINDING_MULTI_HPP
#define BINDING_MULTI_HPP

#include <boost/python.hpp>
#include "Defs.hpp"
#include "peer/converters/tuple.hpp"

namespace binding {

template<typename Timpl>
inline boost::python::objects::class_base register_multi_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    return class_<impl_type, bases<typename impl_type::base_type>,
           boost::shared_ptr<impl_type>, boost::noncopyable>(name,
        init<typename impl_type::identifier_type,
             typename impl_type::simulator_type&,
             Real>())
        .add_property("last_reaction",
            make_function(&impl_type::last_reaction,
                return_value_policy<return_by_value>()))
        .add_property("last_event",
            make_function(&impl_type::last_event,
                return_value_policy<return_by_value>()))
        .add_property("shells",
            make_function(&impl_type::get_shells))
        .def("add_particle", &impl_type::add_particle)
        .def("add_shell", &impl_type::add_shell)
        .def("step", &impl_type::step)
        ;
}

} // namespace binding

#endif /* BINDING_MULTI_HPP */

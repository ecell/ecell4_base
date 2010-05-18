#ifndef BINDING_CYLINDRICAL_SINGLE_HPP
#define BINDING_CYLINDRICAL_SINGLE_HPP

#include <boost/python.hpp>
#include "peer/converters/tuple.hpp"

namespace binding {

template<typename Timpl>
inline boost::python::objects::class_base register_cylindrical_single_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    return class_<impl_type, bases<typename impl_type::base_type>,
           boost::shared_ptr<impl_type> >(name,
        init<typename impl_type::surface_id_type,
             typename impl_type::particle_id_pair,
             typename impl_type::shell_id_pair,
             typename impl_type::reaction_rule_vector const&>())
        .add_property("particle",
            make_function(&impl_type::particle,
                return_value_policy<return_by_value>()))
        .add_property("shell",
            make_function(&impl_type::shell,
                return_value_policy<return_by_value>()))
        .add_property("reactions",
            make_function(&impl_type::reactions,
                return_value_policy<return_by_value>()))
        .add_property("k_tot",
            make_function(&impl_type::k_tot,
                return_value_policy<return_by_value>()))
        ;
}

} // namespace binding

#endif /* BINDING_CYLINDRICAL_SINGLE_HPP */

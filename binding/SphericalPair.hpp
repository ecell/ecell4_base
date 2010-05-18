#ifndef BINDING_SPHERICAL_PAIR_HPP
#define BINDING_SPHERICAL_PAIR_HPP

#include <boost/python.hpp>

namespace binding {

template<typename Timpl>
inline boost::python::objects::class_base register_spherical_pair_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    return class_<impl_type, bases<typename impl_type::base_type>,
           boost::shared_ptr<impl_type> >(name,
        init<typename impl_type::surface_id_type,
             typename impl_type::particle_id_pair,
             typename impl_type::particle_id_pair,
             typename impl_type::shell_id_pair,
             typename impl_type::length_type,
             typename impl_type::length_type>())
        .add_property("shell",
            make_function(&impl_type::shell,
                return_value_policy<return_by_value>()))
        .add_property("r0",
            make_function(&impl_type::r0,
                return_value_policy<return_by_value>()))
        .add_property("rt",
            make_function(&impl_type::rt,
                return_value_policy<return_by_value>()))
        .add_property("a_R",
            make_function(&impl_type::a_R,
                return_value_policy<return_by_value>()))
        .add_property("a_r",
            make_function(&impl_type::a_r,
                return_value_policy<return_by_value>()))
        .add_property("sigma",
            make_function(&impl_type::sigma,
                return_value_policy<return_by_value>()))
        .add_property("D_tot",
            make_function(&impl_type::D_tot,
                return_value_policy<return_by_value>()))
        .add_property("D_geom",
            make_function(&impl_type::D_geom,
                return_value_policy<return_by_value>()))
        .add_property("D_R",
            make_function(&impl_type::D_R,
                return_value_policy<return_by_value>()));
}

} // namespace binding

#endif /* BINDING_SPHERICAL_PAIR_HPP */

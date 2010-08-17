#ifndef BINDING_CYLINDER_HPP
#define BINDING_CYLINDER_HPP

#include <boost/python.hpp>
#include "peer/utils.hpp"

namespace binding {

template<typename Timpl_>
static std::string Cylinder___str__(Timpl_* impl)
{
    return boost::lexical_cast<std::string>(*impl);
}

template<typename Timpl_>
inline boost::python::objects::class_base register_cylinder_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl_ impl_type;

    return class_<impl_type>(name)
        .def(init<typename impl_type::position_type, 
                  typename impl_type::length_type,
                  typename impl_type::position_type, 
                  typename impl_type::length_type>())
        .add_property("position",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::position_type,
                    &impl_type::position,
                    &impl_type::position>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::position_type,
                    &impl_type::position,
                    &impl_type::position>::set))
        .add_property("radius",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::length_type,
                    &impl_type::radius,
                    &impl_type::radius>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::length_type,
                    &impl_type::radius,
                    &impl_type::radius>::set))
        .add_property("size",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::length_type,
                    &impl_type::size,
                    &impl_type::size>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::length_type,
                    &impl_type::size,
                    &impl_type::size>::set))
        .add_property("unit_z",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::position_type,
                    &impl_type::unit_z,
                    &impl_type::unit_z>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::position_type,
                    &impl_type::unit_z,
                    &impl_type::unit_z>::set))
        .def("__str__", &Cylinder___str__<impl_type>)
        .def("show", &impl_type::show);

}

} // namespace binding

#endif /* BINDING_CYLINDER_HPP */

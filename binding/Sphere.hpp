#ifndef BINDING_SPHERE_HPP
#define BINDING_SPHERE_HPP

#include <boost/python.hpp>
#include "peer/utils.hpp"

namespace binding {

template<typename Timpl_>
static boost::python::object Sphere___getitem__(Timpl_ const& obj, int index)
{
    switch (index)
    {
    default:
        PyErr_SetString(PyExc_IndexError, "index out of range");
        boost::python::throw_error_already_set();
        break;
    case 0:
        return boost::python::object(obj.position());
    case 1:
        return boost::python::object(obj.radius());
    }

    return boost::python::object();
}

template<typename Timpl_>
static std::string Sphere___str__(Timpl_* impl)
{
    return boost::lexical_cast<std::string>(*impl);
}

template<typename Timpl_>
inline boost::python::objects::class_base register_sphere_class(char const *name)
{
    using namespace boost::python;
    typedef Timpl_ impl_type;

    return class_<impl_type>(name)
        .def(init<typename impl_type::position_type, 
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
        .def("__getitem__", &Sphere___getitem__<impl_type>)
        .def("__str__", &Sphere___str__<impl_type>)
        .def("show", &impl_type::show);
}

} // namespace binding

#endif /* BINDING_SPHERE_HPP */

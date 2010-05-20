#ifndef BINDING_STRUCTURE_HPP
#define BINDING_STRUCTURE_HPP

#include <boost/python.hpp>

namespace binding {

template<typename Timpl>
inline boost::python::objects::class_base register_structure_class(char const *name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    return class_<impl_type, boost::shared_ptr<impl_type>,
                  boost::noncopyable>(name, no_init)
        .add_property("id", 
            make_function(&impl_type::id,
                          return_value_policy<return_by_value>()))
        .def("random_position", &impl_type::random_position)
        .def("random_vector", &impl_type::random_vector)
        .def("bd_displacement", &impl_type::bd_displacement)
        ;
}

} // namespace binding

#endif /* BINDING_STRUCTURE_HPP */

#ifndef BINDING_STRUCTURE_HPP
#define BINDING_STRUCTURE_HPP

#include <boost/python.hpp>

namespace binding {

template<typename Timpl>
inline boost::python::objects::class_base register_structure_class(char const *name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    return class_<impl_type>(name, no_init)
        .add_property("id", 
            make_function(&impl_type::id,
                          return_value_policy<return_by_value>()))
        ;
}

} // namespace binding

#endif /* BINDING_STRUCTURE_HPP */

#ifndef BINDING_SINGLE_HPP
#define BINDING_SINGLE_HPP

#include <boost/python.hpp>

namespace binding {

template<typename Timpl>
inline boost::python::objects::class_base register_single_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    return class_<impl_type, bases<typename impl_type::base_type>,
           boost::shared_ptr<impl_type> >(name, no_init)
        .add_property("particle",
            make_function(&impl_type::particle,
                return_value_policy<return_by_value>()))
        ;
}

} // namespace binding

#endif /* BINDING_SINGLE_HPP */

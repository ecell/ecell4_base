#ifndef BINDING_SHAPED_DOMAIN_HPP
#define BINDING_SHAPED_DOMAIN_HPP

#include <boost/python.hpp>
#include "peer/utils.hpp"

namespace binding {

template<typename Timpl>
inline boost::python::objects::class_base register_shaped_domain_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    return class_<impl_type, boost::shared_ptr<impl_type>,
                  boost::noncopyable>(name, no_init)
        .add_property("position", 
            make_function(&impl_type::position,
                return_value_policy<return_by_value>()))
        .add_property("size", 
            make_function(&impl_type::size,
                return_value_policy<return_by_value>()))
        ;
}

} // namespace binding

#endif /* BINDING_SHAPED_DOMAIN_HPP */


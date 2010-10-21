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
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::position_type,
                    &impl_type::position, &impl_type::position>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::position_type,
                    &impl_type::position, &impl_type::position>::set,
                return_value_policy<return_by_value>()))
        .add_property("size", 
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::length_type,
                    &impl_type::size, &impl_type::size>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::length_type,
                    &impl_type::size, &impl_type::size>::set,
                return_value_policy<return_by_value>()))
        ;
}

} // namespace binding

#endif /* BINDING_SHAPED_DOMAIN_HPP */


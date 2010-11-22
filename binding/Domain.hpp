#ifndef BINDING_DOMAIN_HPP
#define BINDING_DOMAIN_HPP

#include <boost/python.hpp>
#include "peer/utils.hpp"

namespace binding {

template<typename Timpl>
inline boost::python::objects::class_base register_domain_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    return class_<impl_type, boost::shared_ptr<impl_type>,
                  boost::noncopyable>(name, no_init)
        .add_property("id", 
            make_function(&impl_type::id,
                return_value_policy<return_by_value>()))
        .add_property("last_time", 
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type, typename impl_type::time_type,
                    &impl_type::last_time, &impl_type::last_time>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type, typename impl_type::time_type,
                    &impl_type::last_time, &impl_type::last_time>::set))
        .add_property("event", 
             make_function(
                 &peer::util::reference_accessor_wrapper<
                    impl_type, typename impl_type::event_id_pair_type,
                    &impl_type::event, &impl_type::event>::get,
                 return_value_policy<return_by_value>()),
             make_function(
                 &peer::util::reference_accessor_wrapper<
                    impl_type, typename impl_type::event_id_pair_type,
                    &impl_type::event, &impl_type::event>::set))
        .add_property("__repr__", make_function(&impl_type::as_string))
        ;
}

} // namespace binding

#endif /* BINDING_DOMAIN_HPP */

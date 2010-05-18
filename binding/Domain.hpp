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

    return class_<impl_type>(name, no_init)
        .add_property("surface_id", 
            make_function(&impl_type::surface_id,
                return_value_policy<return_by_value>()))
        .add_property("event_id", 
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type, typename impl_type::event_id_type,
                    &impl_type::event_id, &impl_type::event_id>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type, typename impl_type::event_id_type,
                    &impl_type::event_id, &impl_type::event_id>::set))
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
        .add_property("event_kind", 
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type, typename impl_type::event_kind_type,
                    &impl_type::event_kind, &impl_type::event_kind>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type, typename impl_type::event_kind_type,
                    &impl_type::event_kind, &impl_type::event_kind>::set))
        ;
}

} // namespace binding

#endif /* BINDING_DOMAIN_HPP */

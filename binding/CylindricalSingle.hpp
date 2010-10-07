#ifndef BINDING_CYLINDRICAL_SINGLE_HPP
#define BINDING_CYLINDRICAL_SINGLE_HPP

#include <boost/python.hpp>
#include "peer/converters/tuple.hpp"
#include "fcpair_reference_accessor_wrapper.hpp"

namespace binding {

template<typename Timpl>
inline boost::python::objects::class_base register_cylindrical_single_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    return class_<impl_type, bases<typename impl_type::base_type>,
           boost::shared_ptr<impl_type>, boost::noncopyable>(name,
        init<typename impl_type::identifier_type,
             typename impl_type::particle_id_pair,
             typename impl_type::shell_id_pair>())
        .add_property("shell",
            make_function(
                 &fcpair_reference_accessor_wrapper<
                    impl_type, typename impl_type::shell_id_pair,
                    &impl_type::shell, &impl_type::shell>::get,
                 return_value_policy<return_by_value>()),
            make_function(
                 &fcpair_reference_accessor_wrapper<
                    impl_type, typename impl_type::shell_id_pair,
                    &impl_type::shell, &impl_type::shell>::set))
        ;
}

} // namespace binding

#endif /* BINDING_CYLINDRICAL_SINGLE_HPP */

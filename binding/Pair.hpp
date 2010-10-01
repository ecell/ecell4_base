#ifndef BINDING_PAIR_HPP
#define BINDING_PAIR_HPP

#include <boost/python.hpp>
#include "peer/converters/sequence.hpp"

namespace binding {

template<typename Timpl>
inline boost::python::objects::class_base register_pair_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    peer::converters::register_range_to_tuple_converter<
            typename impl_type::particle_array_type>();

    return class_<impl_type, bases<typename impl_type::base_type>,
           boost::shared_ptr<impl_type>, boost::noncopyable>(name, no_init)
        .add_property("particles",
            make_function(static_cast<typename impl_type::particle_array_type const&(impl_type::*)() const>(&impl_type::particles),
                return_value_policy<return_by_value>()));
}

} // namespace binding

#endif /* BINDING_PAIR_HPP */


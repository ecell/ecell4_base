#ifndef BINDING_PARTICLE_SIMULATOR_HPP
#define BINDING_PARTICLE_SIMULATOR_HPP

#include <boost/python.hpp>
#include <boost/python.hpp>
#include "peer/util/to_native_converter.hpp"
#include "peer/wrappers/generator/generator_wrapper.hpp"
#include "peer/converters/tuple.hpp"

namespace binding {

template<typename Timpl>
void register_particle_simulator_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    class_<impl_type, boost::noncopyable>(name, no_init)
        .add_property("world",
            make_function(
                &impl_type::world,
                return_value_policy<return_by_value>()))
        .add_property("network_rules",
            make_function(
                &impl_type::network_rules,
                return_value_policy<return_by_value>()))
        .add_property("reaction_recorder",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    boost::shared_ptr<typename impl_type::reaction_recorder_type>,
                    &impl_type::reaction_recorder, &impl_type::reaction_recorder>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    boost::shared_ptr<typename impl_type::reaction_recorder_type>,
                    &impl_type::reaction_recorder, &impl_type::reaction_recorder>::set))
        .add_property("rng",
            make_function(
                &impl_type::rng,
                return_value_policy<reference_existing_object>()))
        .add_property("paranoiac",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type, bool,
                    &impl_type::paranoiac, &impl_type::paranoiac>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type, bool,
                    &impl_type::paranoiac, &impl_type::paranoiac>::set))
        .add_property("t", &impl_type::t)
        .add_property("dt", &impl_type::dt)
        .add_property("num_steps", &impl_type::num_steps)
        .def("step", static_cast<void(impl_type::*)()>(&impl_type::step))
        .def("step", static_cast<bool(impl_type::*)(typename impl_type::time_type)>(&impl_type::step))
        ;
}

} // namespace binding
#endif /* BINDING_PARTICLE_SIMULATOR_HPP */

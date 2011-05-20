#ifndef BINDING_BD_PROPAGATOR_HPP
#define BINDING_BD_PROPAGATOR_HPP

#include <boost/python.hpp>
#include "utils/pair.hpp"
#include "peer/utils.hpp"
#include "peer/wrappers/range/pyiterable_range.hpp"
#include "peer/converters/sequence.hpp"

namespace binding {

template<typename Timpl_>
static void BDPropagator_propagate_all(Timpl_& self)
{
    while (self());
}

template<typename Timpl_>
inline boost::python::objects::class_base register_bd_propagator_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl_ impl_type;
    typedef typename impl_type::traits_type simulator_traits_type;
    typedef typename simulator_traits_type::world_type world_type;
    typedef typename world_type::particle_id_type particle_id_type;
    typedef typename simulator_traits_type::network_rules_type network_rules_type;
    typedef typename world_type::traits_type::rng_type rng_type;
    typedef typename simulator_traits_type::time_type time_type;
    typedef typename world_type::particle_container_type particle_container_type;
    typedef typename impl_type::reaction_recorder_type reaction_recorder_type;
    typedef typename impl_type::volume_clearer_type volume_clearer_type;

    peer::converters::register_pyiterable_range_converter<particle_id_type>();
    return class_<impl_type, boost::noncopyable>(
        name, init<
            particle_container_type&, network_rules_type const&, rng_type&,
            time_type, int, reaction_recorder_type*, volume_clearer_type*,
            peer::wrappers::pyiterable_range<particle_id_type> >())
        .def(init<
            particle_container_type&, network_rules_type const&, rng_type&,
            time_type, int, reaction_recorder_type*, volume_clearer_type*,
            typename get_select_first_range<typename world_type::particle_id_pair_range>::type>())
        .add_property("rejected_move_count",
            &impl_type::get_rejected_move_count)
        .def("__call__", &impl_type::operator())
        .def("propagate_all", &BDPropagator_propagate_all<impl_type>)
        ;
}

} // namespace binding

#endif /* BINDING_BD_PROPAGATOR_HPP */

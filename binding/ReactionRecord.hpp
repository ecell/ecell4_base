#ifndef BINDING_REACTION_RECORD_HPP
#define BINDING_REACTION_RECORD_HPP

#include <boost/python.hpp>
#include "peer/converters/sequence.hpp"

namespace binding {

template<typename Timpl_>
boost::python::objects::class_base register_reaction_record_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl_ impl_type;

    peer::converters::register_range_to_tuple_converter<typename impl_type::reactants_type>();
    peer::converters::register_range_to_tuple_converter<typename impl_type::products_type>();
    peer::converters::register_pyiterable_range_converter<typename impl_type::particle_id_type>();

    return class_<impl_type>(name, init<>())
        .def(init<typename impl_type::reaction_rule_id_type const&,
                  peer::wrappers::pyiterable_range<
                    typename impl_type::particle_id_type>,
                  typename impl_type::particle_id_type const&>())
        .def(init<typename impl_type::reaction_rule_id_type const&,
                  peer::wrappers::pyiterable_range<
                    typename impl_type::particle_id_type>,
                  typename impl_type::particle_id_type const&,
                  typename impl_type::particle_id_type const&>())
        .def(self == self)
        .def(self != self)
        .add_property("reaction_rule_id",
            make_function(&impl_type::reaction_rule_id,
                          return_value_policy<return_by_value>()))
        .add_property("reactants",
            make_function(&impl_type::reactants,
                          return_value_policy<return_by_value>()))
        .add_property("products",
            make_function(&impl_type::products,
                          return_value_policy<return_by_value>()));
}

} // namespace binding

#endif /* BINDING_REACTION_RECORD_HPP */

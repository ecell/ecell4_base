#ifndef BINDING_REACTION_RULE_INFO_HPP
#define BINDING_REACTION_RULE_INFO_HPP

#include <boost/python.hpp>
#include "peer/converters/sequence.hpp"

namespace binding {

template<typename Timpl>
inline void register_reaction_rule_info_class(char const *name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    class_<impl_type>(name,
        init<typename impl_type::identifier_type,
             typename impl_type::rate_type,
             twofold_container<typename impl_type::species_id_type>,
             std::vector<typename impl_type::species_id_type> >())
        .add_property("id", 
            make_function(&impl_type::id,
                          return_value_policy<return_by_value>()))
        .add_property("k", make_function(&impl_type::k))
        .add_property("products",
            make_function(&impl_type::get_products,
                          return_value_policy<return_by_value>()))
        .add_property("reactants",
            make_function(&impl_type::get_reactants,
                          return_value_policy<return_by_value>()));

    peer::converters::register_range_to_tuple_converter<typename impl_type::species_id_range>();

    peer::converters::register_iterable_to_range_converter<typename impl_type::species_id_range>();
}

} // namespace binding

#endif /* BINDING_REACTION_RULE_INFO_HPP */

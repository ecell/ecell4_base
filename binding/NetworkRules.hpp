#ifndef BINDING_NETWORK_RULES_HPP
#define BINDING_NETWORK_RULES_HPP
  
#include <boost/python.hpp>

namespace binding {

template<typename Timpl_>
boost::python::objects::class_base register_network_rules_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl_ impl_type;

    return class_<impl_type, boost::noncopyable>(name, no_init)
        .def("add_reaction_rule", &impl_type::add_reaction_rule)
        .def("remove_reaction_rule", &impl_type::remove_reaction_rule)
        .def("query_reaction_rule", static_cast<typename impl_type::reaction_rule_generator*(impl_type::*)(typename impl_type::species_id_type const&) const>(&impl_type::query_reaction_rule), return_value_policy<return_by_value>())
        .def("query_reaction_rule", static_cast<typename impl_type::reaction_rule_generator*(impl_type::*)(typename impl_type::species_id_type const&, typename impl_type::species_id_type const&) const>(&impl_type::query_reaction_rule), return_value_policy<return_by_value>())
        ;
}

} // namespace binding

#endif /* BINDING_NETWORK_RULES_HPP */

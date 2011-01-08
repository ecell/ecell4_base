#ifndef BINDING_NETWORK_RULES_WRAPPER_HPP
#define BINDING_NETWORK_RULES_WRAPPER_HPP

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/python.hpp>
#include "../NetworkRules.hpp"

namespace binding {

template<typename Timpl>
boost::python::objects::class_base register_network_rules_wrapper_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    return class_<impl_type, boost::shared_ptr<impl_type>,
                  boost::noncopyable>(name, init<NetworkRules const&>())
        .def("query_reaction_rule", (typename impl_type::reaction_rule_vector const&(impl_type::*)(typename impl_type::species_id_type const&) const)&impl_type::query_reaction_rule,
            return_value_policy<return_by_value>())
        .def("query_reaction_rule", (typename impl_type::reaction_rule_vector const&(impl_type::*)(typename impl_type::species_id_type const&, typename impl_type::species_id_type const&) const)&impl_type::query_reaction_rule,
            return_value_policy<return_by_value>())
        ;
}

} // namespace binding

#endif /* BINDING_NETWORK_RULES_WRAPPER_HPP */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/lexical_cast.hpp>
#include "utils.hpp"
#include "exceptions.hpp"
#include "basic_network_rules_impl.hpp"

basic_network_rules_impl::~basic_network_rules_impl()
{
}

basic_network_rules_impl::basic_network_rules_impl()
{
}

void basic_network_rules_impl::add_reaction_rule(reaction_rule const& r)
{
    if (!reaction_rules_map_[r.get_reactants()].insert(r).second)
        throw already_exists(boost::lexical_cast<std::string>(r));
}
    
basic_network_rules_impl::reaction_rule_generator*
basic_network_rules_impl::query_reaction_rule(species_type const* r1)
{
    reaction_rules_map::const_iterator i(
            reaction_rules_map_.find(reaction_rule::reactants(r1)));
    if (i == reaction_rules_map_.end())
    {
        return 0;
    }
    return make_range_generator((*i).second);
}

basic_network_rules_impl::reaction_rule_generator*
basic_network_rules_impl::query_reaction_rule(species_type const* r1, species_type const* r2)
{
    reaction_rules_map::const_iterator i(
            reaction_rules_map_.find(reaction_rule::reactants(r1, r2)));
    if (i == reaction_rules_map_.end())
    {
        return 0;
    }
    return make_range_generator((*i).second);
}

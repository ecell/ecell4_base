#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/lexical_cast.hpp>
#include "utils.hpp"
#include "exceptions.hpp"
#include "BasicNetworkRulesImpl.hpp"

BasicNetworkRulesImpl::~BasicNetworkRulesImpl()
{
}

BasicNetworkRulesImpl::BasicNetworkRulesImpl()
{
}

void BasicNetworkRulesImpl::add_reaction_rule(ReactionRule const& r)
{
    if (!reaction_rules_map_[r.get_reactants()].insert(r).second)
        throw already_exists(boost::lexical_cast<std::string>(r));
}
    
BasicNetworkRulesImpl::reaction_rule_generator*
BasicNetworkRulesImpl::query_reaction_rule(SpeciesType const* r1)
{
    reaction_rules_map::const_iterator i(
            reaction_rules_map_.find(ReactionRule::Reactants(r1)));
    if (i == reaction_rules_map_.end())
    {
        return 0;
    }
    return make_range_generator((*i).second);
}

BasicNetworkRulesImpl::reaction_rule_generator*
BasicNetworkRulesImpl::query_reaction_rule(SpeciesType const* r1, SpeciesType const* r2)
{
    reaction_rules_map::const_iterator i(
            reaction_rules_map_.find(ReactionRule::Reactants(r1, r2)));
    if (i == reaction_rules_map_.end())
    {
        return 0;
    }
    return make_range_generator((*i).second);
}

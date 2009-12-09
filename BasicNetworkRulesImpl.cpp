#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/lexical_cast.hpp>
#include "utils/range_support.hpp"
#include "exceptions.hpp"
#include "generator.hpp"
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
BasicNetworkRulesImpl::query_reaction_rule(SpeciesTypeID const& r1) const
{
    reaction_rules_map::const_iterator i(
            reaction_rules_map_.find(ReactionRule::Reactants(r1)));
    if (i == reaction_rules_map_.end())
    {
        return 0;
    }
    return make_range_generator<ReactionRule>((*i).second);
}

BasicNetworkRulesImpl::reaction_rule_generator*
BasicNetworkRulesImpl::query_reaction_rule(SpeciesTypeID const& r1, SpeciesTypeID const& r2) const
{
    reaction_rules_map::const_iterator i(
            reaction_rules_map_.find(ReactionRule::Reactants(r1, r2)));
    if (i == reaction_rules_map_.end())
    {
        return 0;
    }
    return make_range_generator<ReactionRule>((*i).second);
}

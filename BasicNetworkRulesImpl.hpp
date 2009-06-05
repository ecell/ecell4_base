#ifndef BASIC_NETWORK_RULES_IMPL_HPP
#define BASIC_NETWORK_RULES_IMPL_HPP

#include <map>
#include <set>

#include "NetworkRules.hpp"
#include "ReactionRule.hpp"

class BasicNetworkRulesImpl: public NetworkRules
{
    typedef std::set<ReactionRule> reaction_rule_set;
    typedef std::map<ReactionRule::Reactants, reaction_rule_set> reaction_rules_map;

public:
    virtual void add_reaction_rule(ReactionRule const&);

    virtual reaction_rule_generator* query_reaction_rule(SpeciesType const* r1) const;

    virtual reaction_rule_generator* query_reaction_rule(SpeciesType const* r1, SpeciesType const* r2) const;

    virtual ~BasicNetworkRulesImpl();

    BasicNetworkRulesImpl();

private:
    reaction_rules_map reaction_rules_map_;
};

#endif /* BASIC_NETWORK_RULES_IMPL_HPP */

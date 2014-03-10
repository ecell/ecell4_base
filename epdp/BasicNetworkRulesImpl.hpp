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
    typedef ReactionRule::identifier_type identifier_type;

public:
    virtual identifier_type add_reaction_rule(ReactionRule const&);
    virtual void remove_reaction_rule(ReactionRule const&);

    virtual reaction_rule_generator* query_reaction_rule(SpeciesTypeID const& r1) const;

    virtual reaction_rule_generator* query_reaction_rule(SpeciesTypeID const& r1, SpeciesTypeID const& r2) const;

    virtual ~BasicNetworkRulesImpl();

    BasicNetworkRulesImpl();

private:
    reaction_rules_map reaction_rules_map_;
    identifier_type serial_;
};

#endif /* BASIC_NETWORK_RULES_IMPL_HPP */

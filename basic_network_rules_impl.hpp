#ifndef BASIC_NETWORK_RULES_IMPL_HPP
#define BASIC_NETWORK_RULES_IMPL_HPP

#include <map>
#include <set>

#include "network_rules.hpp"
#include "reaction_rule.hpp"

class basic_network_rules_impl: public network_rules
{
    typedef std::set<reaction_rule> reaction_rule_set;
    typedef std::map<reaction_rule::reactants, reaction_rule_set> reaction_rules_map;

public:
    virtual void add_reaction_rule(reaction_rule const&);

    virtual reaction_rule_generator* query_reaction_rule(species_type const* r1);

    virtual reaction_rule_generator* query_reaction_rule(species_type const* r1, species_type const* r2);

    virtual ~basic_network_rules_impl();

    basic_network_rules_impl();

private:
    reaction_rules_map reaction_rules_map_;
};

#endif /* BASIC_NETWORK_RULES_IMPL_HPP */

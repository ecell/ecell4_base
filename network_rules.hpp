#ifndef NETWORK_RULES_HPP
#define NETWORK_RULES_HPP

#include <map>
#include <boost/range/iterator_range.hpp>
#include "reaction_rule.hpp"

class network_rules
{
public:
    virtual void add_reaction_rule(reaction_rule const&) = 0;

    virtual reaction_rule_generator* query_reaction_rule(species_type const* r1) = 0;

    virtual reaction_rule_generator* query_reaction_rule(species_type const* r1, species_type const* r2) = 0;

    virtual ~network_rules() = 0;
};

#endif /* NETWORK_RULES_HPP */

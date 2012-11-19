#ifndef __REACTION_RULE_HPP
#define __REACTION_RULE_HPP

#include <vector>

#include "types.hpp"
#include "Species.hpp"


namespace ecell4
{

class ReactionRule
{
public:

    typedef std::vector<Species> SpeciesVector;

    virtual SpeciesVector const& reactants() const = 0;
    virtual SpeciesVector const& products() const = 0;

    virtual void add_product(Species const& sp) const = 0;
};

typedef std::vector<ReactionRule> ReactionRuleVector;

} // ecell4

#endif /* __REACTION_RULE_HPP */

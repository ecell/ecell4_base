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

    ReactionRule()
        : reactants_(), products_()
    {
        ;
    }

    Real k() const
    {
        return k_;
    }

    virtual SpeciesVector const& reactants() const
    {
        return reactants_;
    }

    virtual SpeciesVector const& products() const
    {
        return products_;
    }

    virtual void add_product(Species const& sp) const
    {
        ;
    }

protected:

    Real k_;
    SpeciesVector reactants_, products_;
};

typedef std::vector<ReactionRule> ReactionRuleVector;

} // ecell4

#endif /* __REACTION_RULE_HPP */

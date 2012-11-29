#ifndef __REACTION_RULE_HPP
#define __REACTION_RULE_HPP

#include <vector>
#include <stdexcept>

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

    SpeciesVector const& reactants() const
    {
        return reactants_;
    }

    SpeciesVector const& products() const
    {
        return products_;
    }

    void set_k(Real const& k)
    {
        if (k <= 0)
        {
            throw std::invalid_argument("a kinetic rate must be positive.");
        }
        k_ = k;
    }

    void add_reactant(Species const& sp)
    {
        reactants_.push_back(sp);
    }

    void add_product(Species const& sp)
    {
        products_.push_back(sp);
    }

protected:

    Real k_;
    SpeciesVector reactants_, products_;
};

typedef std::vector<ReactionRule> ReactionRuleVector;

} // ecell4

#endif /* __REACTION_RULE_HPP */

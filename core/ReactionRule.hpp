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

    /**
     * a type of the container of reactants
     * std::multiset allows multiple keys with equal values,
     * but looses the original order at the registration.
     */
    typedef std::multiset<Species> reactants_type;
    typedef std::multiset<Species> products_type;

    ReactionRule()
        : k_(0), reactants_(), products_()
    {
        ;
    }

    Real k() const
    {
        return k_;
    }

    reactants_type const& reactants() const
    {
        return reactants_;
    }

    products_type const& products() const
    {
        return products_;
    }

    void set_k(Real const& k)
    {
        if (k < 0)
        {
            throw std::invalid_argument("a kinetic rate must be positive.");
        }
        k_ = k;
    }

    void add_reactant(Species const& sp)
    {
        reactants_.insert(sp);
    }

    void add_product(Species const& sp)
    {
        products_.insert(sp);
    }

protected:

    Real k_;
    reactants_type reactants_;
    products_type products_;
};

inline bool operator<(ReactionRule const& lhs, ReactionRule const& rhs)
{
    if (lhs.reactants() < rhs.reactants())
    {
        return true;
    }
    else if (lhs.reactants() > rhs.reactants())
    {
        return false;
    }
    return (lhs.products() < rhs.products());
}

inline bool operator==(ReactionRule const& lhs, ReactionRule const& rhs)
{
    return ((lhs.reactants() == rhs.reactants())
            && (lhs.products() == rhs.products()));
}

inline bool operator!=(ReactionRule const& lhs, ReactionRule const& rhs)
{
    return !(lhs == rhs);
}

typedef std::vector<ReactionRule> ReactionRuleVector;

} // ecell4

#endif /* __REACTION_RULE_HPP */

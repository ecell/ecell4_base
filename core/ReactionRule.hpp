#ifndef __REACTION_RULE_HPP
#define __REACTION_RULE_HPP

#include <set>
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
    typedef std::multiset<Species> reactant_container_type;
    typedef std::multiset<Species> product_container_type;

    ReactionRule()
        : k_(0), reactants_(), products_()
    {
        ;
    }

    Real k() const
    {
        return k_;
    }

    reactant_container_type const& reactants() const
    {
        return reactants_;
    }

    product_container_type const& products() const
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
    reactant_container_type reactants_;
    product_container_type products_;
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

} // ecell4

#endif /* __REACTION_RULE_HPP */

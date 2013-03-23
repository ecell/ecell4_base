#ifndef __ECELL4_REACTION_RULE_HPP
#define __ECELL4_REACTION_RULE_HPP

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
     * when changing this type into the ordered one,
     * please modify NetworkModel too.
     */
    typedef std::multiset<Species> reactant_container_type;
    typedef std::multiset<Species> product_container_type;

public:

    ReactionRule()
        : k_(0), reactants_(), products_()
    {
        ;
    }

    Real k() const
    {
        return k_;
    }

    const reactant_container_type& reactants() const
    {
        return reactants_;
    }

    const product_container_type& products() const
    {
        return products_;
    }

    void set_k(const Real& k)
    {
        if (k < 0)
        {
            throw std::invalid_argument("a kinetic rate must be positive.");
        }
        k_ = k;
    }

    void add_reactant(const Species& sp)
    {
        reactants_.insert(sp);
    }

    void add_product(const Species& sp)
    {
        products_.insert(sp);
    }

protected:

    Real k_;
    reactant_container_type reactants_;
    product_container_type products_;
};

inline bool operator<(const ReactionRule& lhs, const ReactionRule& rhs)
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

inline bool operator==(const ReactionRule& lhs, const ReactionRule& rhs)
{
    return ((lhs.reactants() == rhs.reactants())
            && (lhs.products() == rhs.products()));
}

inline bool operator!=(const ReactionRule& lhs, const ReactionRule& rhs)
{
    return !(lhs == rhs);
}

} // ecell4

#endif /* __ECELL4_REACTION_RULE_HPP */

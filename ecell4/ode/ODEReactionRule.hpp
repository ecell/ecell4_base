#ifndef __ECELL4_ODE_REACTION_RULE_HPP
#define __ECELL4_ODE_REACTION_RULE_HPP

// #include <set>
#include <stdexcept>
#include <map>
#include <vector>

#include <ecell4/core/types.hpp>
#include <ecell4/core/Species.hpp>
//#include <ecell4/core/Ratelaw.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include "ODERatelaw.hpp"


namespace ecell4
{

namespace ode
{

class ODEReactionRule
{
public:

    /**
     * a type of the container of reactants
     * std::multiset allows multiple keys with equal values,
     * but looses the original order at the registration.
     * when changing this type into the ordered one,
     * please modify NetworkModel too.
     */
    // To be compatible with ReactionRule of core.
    typedef std::vector<Species> reactant_container_type;
    typedef std::vector<Species> product_container_type;

    typedef std::vector<Real> coefficient_container_type;

    typedef std::vector<std::pair<Real, Species> > reaction_leftside_container_type;
    typedef std::vector<std::pair<Real, Species> > reaction_rightside_container_type;

public:

    ODEReactionRule()
        : k_(0), reactants_(), products_()
    {
        ;
    }
    ODEReactionRule(const ODEReactionRule &rr)
        : k_(rr.k_), reactants_(rr.reactants_), products_(rr.products_)
    {
        if(rr.has_ratelaw())
        {
            this->set_ratelaw(rr.get_ratelaw());
        }
    }

    ODEReactionRule(
        const reaction_leftside_container_type& reactants,
        const reaction_rightside_container_type& products,
        const Real &k = 0.0)
        : k_(k), reactants_(reactants), products_(products)
    {
        ;
    }

    ODEReactionRule(const ecell4::ReactionRule& rr)
        : k_(rr.k())
    {
        // This constructor is for compatibility with ReactionRule defined in core-module.
        for(reactant_container_type::const_iterator it(rr.reactants().begin()); 
                it != rr.reactants().end(); it++)
        {
            this->add_reactant(*it, 1.0);
        }
        for(product_container_type::const_iterator it(rr.products().begin());
                it != rr.products().end(); it++)
        {
            this->add_product(*it, 1.0);
        }
    }

    ODEReactionRule(
        const reactant_container_type& reactants,
        const product_container_type& products,
        const Real& k = 0.0)
        : k_(k)
    {
        // This constructor is for compatibility with ReactionRule defined in core-module.
        for(reactant_container_type::const_iterator it(reactants.begin()); 
                it != reactants.end(); it++)
        {
            this->add_reactant(*it, 1.0);
        }
        for(product_container_type::const_iterator it(products.begin());
                it != products.end(); it++)
        {
            this->add_product(*it, 1.0);
        }
    }

    Real k() const
    {
        return k_;
    }

    const reactant_container_type reactants() const
    {
        reactant_container_type result;
        for(reaction_leftside_container_type::const_iterator it(reactants_.begin() ); 
                it != reactants_.end(); it++)
        {
            result.push_back(it->second);
        }
        return result;
    }
    const coefficient_container_type reactants_coefficients() const
    {
        coefficient_container_type result;
        for(reaction_leftside_container_type::const_iterator it(reactants_.begin());
                it != reactants_.end(); it++)
        {
            result.push_back(it->first);
        }
        return result;
    }

    const product_container_type products() const
    {
        product_container_type result;
        for(reaction_rightside_container_type::const_iterator it(products_.begin() ); 
                it != products_.end(); it++)
        {
            result.push_back(it->second);
        }
        return result;
    }
    const coefficient_container_type products_coefficients() const
    {
        coefficient_container_type result;
        for(reaction_rightside_container_type::const_iterator it(products_.begin());
                it != products_.end(); it++)
        {
            result.push_back(it->first);
        }
        return result;
    }

    void set_k(const Real& k)
    {
        if (k < 0)
        {
            throw std::invalid_argument("a kinetic rate must be positive.");
        }
        k_ = k;
    }

    void add_reactant(const Species& sp, Real coefficient = 1.0)
    {
        this->reactants_.push_back(std::pair<Real, Species>(coefficient, sp));
    }

    void set_reactant_coefficient(const std::size_t num, const Real new_coeff)
    {
        this->reactants_[num].first = new_coeff;
    }


    void add_product(const Species& sp, Real coefficient = 1.0)
    {
        this->products_.push_back(std::pair<Real, Species>(coefficient, sp));
    }

    void set_product_coefficient(const std::size_t num, const Real new_coeff)
    {
        this->products_[num].first = new_coeff;
    }

    const std::string as_string() const;
    //Integer count(const reactant_container_type& reactants) const;
    //std::vector<ReactionRule> generate(const reactant_container_type& reactants) const;

    /** Ratelaw related functions.
      */

    void set_ratelaw(const boost::shared_ptr<ODERatelaw> ratelaw)
    {
        this->ratelaw_ = ratelaw;
    }

    boost::shared_ptr<ODERatelaw> get_ratelaw() const
    {
        return this->ratelaw_.lock();
    }

    bool has_ratelaw() const
    {
        return !(this->ratelaw_.expired());
    }

    

protected:

    Real k_;
    //reactant_container_type reactants_;
    //product_container_type products_;
    reaction_leftside_container_type reactants_;
    reaction_rightside_container_type products_;

    boost::weak_ptr<ODERatelaw> ratelaw_;
};

#if 0
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
#endif

} // ode

} // ecell4

#endif /* __ECELL4_ODE_REACTION_RULE_HPP */

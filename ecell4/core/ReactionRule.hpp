#ifndef ECELL4_REACTION_RULE_HPP
#define ECELL4_REACTION_RULE_HPP

// #include <set>
#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#include "types.hpp"
#include "Species.hpp"
#include "pyhandler.hpp"
//#include "Ratelaw.hpp"


namespace ecell4
{

class ReactionRuleDescriptor
{
public:
    bool is_available() const
    {   return false;    }
private:

};

class ReactionRuleDescriptorPyfunc
    : public ReactionRuleDescriptor 
{

    typedef std::vector<Species> reactant_container_type;
    typedef std::vector<Species> product_container_type;

    typedef ReactionRuleDescriptor base_type;
    typedef void *pyfunc_type;
    typedef bool (*stepladder_type_rrdescriptor)(
            pyfunc_type, reactant_container_type, product_container_type,  boost::shared_ptr<PyObjectHandler> py_handler );

public:
    bool is_available() const 
    {
        if (pyfunc_ != 0) {return true;}
        return false;
    }
    ReactionRuleDescriptorPyfunc(
            stepladder_type_rrdescriptor stepladder, pyfunc_type pyfunc, boost::shared_ptr<PyObjectHandler> py_handler)
        :pyfunc_(pyfunc)
    {;}
    Real flux() const 
    {
        return 0.;
    }
private:
    pyfunc_type pyfunc_;
};

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
    typedef std::vector<Species> reactant_container_type;
    typedef std::vector<Species> product_container_type;

public:

    enum policy_type
    {
        STRICT = 1L << 0,
        IMPLICIT = 1L << 1,
        DESTROY = 1L << 2
    };

public:

    ReactionRule()
        : k_(0), reactants_(), products_(), policy_(STRICT)
    {
        ;
    }

    ReactionRule(
        const reactant_container_type& reactants,
        const product_container_type& products)
        : k_(0), reactants_(reactants), products_(products), policy_(STRICT)
    {
        ;
    }

    ReactionRule(
        const reactant_container_type& reactants,
        const product_container_type& products,
        const Real& k)
        : k_(k), reactants_(reactants), products_(products), policy_(STRICT)
    {
        ;
    }

    ReactionRule(
        const ReactionRule& rr)
        : k_(rr.k()), reactants_(rr.reactants()), products_(rr.products()), policy_(rr.policy())
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
        reactants_.push_back(sp);
    }

    void add_product(const Species& sp)
    {
        products_.push_back(sp);
    }

    const policy_type policy() const
    {
        return policy_;
    }

    void set_policy(const policy_type policy)
    {
        policy_ = policy;
    }

    const std::string as_string() const;

    inline Integer count(const reactant_container_type& reactants) const
    {
        return this->generate(reactants).size();
    }

    std::vector<ReactionRule> generate(const reactant_container_type& reactants) const;

    /** Ratelaw related functions.
      */
    /*
    void set_ratelaw(const boost::shared_ptr<Ratelaw> ratelaw)
    {
        this->ratelaw_ = ratelaw;
    }

    boost::shared_ptr<Ratelaw> get_ratelaw() const
    {
        return this->ratelaw_.lock();
    }

    bool has_ratelaw() const
    {
        return !(this->ratelaw_.expired());
    }*/
    
    /** ReactionRule Descriptor related functions.
      */
    void set_descriptor(const boost::shared_ptr<ReactionRuleDescriptor> rrd) 
    {
        this->rr_descriptor_ = rrd;
    }
    bool has_descriptor() const 
    {
        return !(this->rr_descriptor_.expired());
    }
    boost::shared_ptr<ReactionRuleDescriptor> get_descriptor() const
    {
        return this->rr_descriptor_.lock();
    }

protected:

    Real k_;
    reactant_container_type reactants_;
    product_container_type products_;

    policy_type policy_;
    //boost::weak_ptr<Ratelaw> ratelaw_;
    boost::weak_ptr<ReactionRuleDescriptor> rr_descriptor_;
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

ReactionRule format_reaction_rule_with_nosort(const ReactionRule& rr);
ReactionRule format_reaction_rule(const ReactionRule& rr);

} // ecell4

#endif /* ECELL4_REACTION_RULE_HPP */

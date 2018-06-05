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
//
#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_io.hpp"


namespace ecell4
{

class ReactionRuleDescriptor
{
public:
    typedef std::vector<Species> reactant_container_type;
    typedef std::vector<Species> product_container_type;
    typedef std::vector<Real> reaction_coefficient_list_type;

    typedef boost::tuple<bool, int, int> descriptor_attribute;  // has_propensity_func, number of reactant coefficients, number of product_coefficients;
public:
    virtual
    bool is_available() const
    {   return true;    }

    descriptor_attribute attribute(void) const
    {   return boost::make_tuple(this->is_available(), reactant_coefficients_.size(), product_coefficients_.size());    }

    virtual
    Real propensity(const std::vector<Real> &r, const std::vector<Real> &p, Real time) const = 0;

    // Accessor of coefficients;
    const reaction_coefficient_list_type &reactant_coefficients(void) const
    {   return this->reactant_coefficients_;    }

    const reaction_coefficient_list_type &product_coefficients(void) const
    {   return this->product_coefficients_; }

    void set_reactant_coefficients(const reaction_coefficient_list_type &new_reactant_coefficients) 
    {
        this->reactant_coefficients_.clear();
        for(int i = 0; i < new_reactant_coefficients.size(); i++) {
            this->reactant_coefficients_.push_back(new_reactant_coefficients[i]);
        }
    }

    void set_product_coefficients(const reaction_coefficient_list_type &new_product_coefficients) 
    {
        this->product_coefficients_.clear();
        for(int i = 0; i < new_product_coefficients.size(); i++) {
            this->product_coefficients_.push_back(new_product_coefficients[i]);
        }
    }

    bool has_coefficients(void) const
    {
        return !(this->reactant_coefficients_.empty() && this->product_coefficients_.empty());
    }
    
private:
    reaction_coefficient_list_type reactant_coefficients_;
    reaction_coefficient_list_type product_coefficients_;
};

class ReactionRuleDescriptorCPPfunc
    : public ReactionRuleDescriptor
{
public:
    typedef Real (*func_type)(const std::vector<Real> &r, const std::vector<Real> &p, Real t);

    ReactionRuleDescriptorCPPfunc(func_type pf) : pf_(pf) {;}

    bool is_available() const
    {
        if (this->pf_ != 0) {return true;}
        return false;
    }

    virtual
    Real propensity(const std::vector<Real> &reactants, const std::vector<Real> &products, Real time) const
    {
        if (this->is_available() ) {
            return (this->pf_)(reactants, products, time);
        } else {
            throw IllegalState("Pointer to the user-defined propensity function is NULL");
        }
    }

private:
    func_type pf_;
};

class ReactionRuleDescriptorPyfunc
    : public ReactionRuleDescriptor 
{
public:
    typedef ReactionRuleDescriptor base_type;
    typedef void *pyfunc_type;
    typedef Real (*stepladder_type_rrdescriptor)(
            pyfunc_type, std::vector<Real>, std::vector<Real>, Real t);

public:
    bool is_available() const 
    {
        if (pyfunc_ != 0 && stepladder_ != 0) {
            return true;
        } else {
            return false;
        }
    }
    ReactionRuleDescriptorPyfunc(
            stepladder_type_rrdescriptor stepladder, pyfunc_type pyfunc, boost::shared_ptr<PyObjectHandler> py_handler)
        :stepladder_(stepladder), pyfunc_(pyfunc), pyobject_handler_(py_handler)
    {
        this->pyobject_handler_->inc_ref(this->pyfunc_);
    }
    virtual
    ~ReactionRuleDescriptorPyfunc()
    {
        this->pyobject_handler_->dec_ref(this->pyfunc_);
    }

    virtual
    Real propensity(const std::vector<Real> &r, const std::vector<Real> &p, Real time) const
    {
        if (stepladder_ == NULL) {  throw IllegalState("stepladder is not registered"); }
        if (pyfunc_ == NULL) {  throw IllegalState("pyfunc is not registered"); }
        Real ret = stepladder_(this->pyfunc_, r, p, time);
        return ret;
    }
private:
    pyfunc_type pyfunc_;
    stepladder_type_rrdescriptor stepladder_;
    boost::shared_ptr<PyObjectHandler> pyobject_handler_;
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
    Real propensity(const std::vector<Real> &r, const std::vector<Real> &p, Real t) const
    {
        if (!has_descriptor()) {
            throw IllegalState("ReactionRule Descriptor has not been registered");
        }
        return this->get_descriptor()->propensity(r, p, t);
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

ReactionRule create_unimolecular_reaction_rule(
    const Species& reactant1, const Species& product1, const Real& k);

ReactionRule create_binding_reaction_rule(
    const Species& reactant1, const Species& reactant2, const Species& product1,
    const Real& k);

ReactionRule create_unbinding_reaction_rule(
    const Species& reactant1, const Species& product1, const Species& product2,
    const Real& k);

ReactionRule create_degradation_reaction_rule(
    const Species& reactant1, const Real& k);

ReactionRule create_synthesis_reaction_rule(
    const Species& product1, const Real& k);

// ReactionRule create_repulsive_reaction_rule(
//     const Species& reactant1, const Species& reactant2);


} // ecell4

#endif /* ECELL4_REACTION_RULE_HPP */

#ifndef ECELL4_REACTION_RULE_DESCRIPTOR_HPP
#define ECELL4_REACTION_RULE_DESCRIPTOR_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#include "types.hpp"
#include "Species.hpp"
#include "pyhandler.hpp"

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

    // The following stores propensity_func, number of reactant coefficients, number of product_coefficients
    typedef boost::tuple<bool, int, int> descriptor_attribute;

public:

    virtual bool is_available() const
    {
        return true;
    }

    descriptor_attribute attribute(void) const
    {
        return boost::make_tuple(
            this->is_available(), reactant_coefficients_.size(), product_coefficients_.size());
    }

    virtual Real propensity(const std::vector<Real> &r, const std::vector<Real> &p, Real time) const = 0;

    // Accessor of coefficients;
    const reaction_coefficient_list_type &reactant_coefficients(void) const
    {
        return this->reactant_coefficients_;
    }

    const reaction_coefficient_list_type &product_coefficients(void) const
    {
        return this->product_coefficients_;
    }

    void set_reactant_coefficients(const reaction_coefficient_list_type &new_reactant_coefficients)
    {
        this->reactant_coefficients_.clear();
        for(int i = 0; i < new_reactant_coefficients.size(); i++)
        {
            this->reactant_coefficients_.push_back(new_reactant_coefficients[i]);
        }
    }

    void set_product_coefficients(const reaction_coefficient_list_type &new_product_coefficients)
    {
        this->product_coefficients_.clear();
        for(int i = 0; i < new_product_coefficients.size(); i++)
        {
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

    ReactionRuleDescriptorCPPfunc(func_type pf)
        : pf_(pf)
    {
        ;
    }

    bool is_available() const
    {
        return (this->pf_ != 0);
    }

    virtual Real propensity(const std::vector<Real> &reactants, const std::vector<Real> &products, Real time) const
    {
        if (this->is_available())
        {
            return (this->pf_)(reactants, products, time);
        }
        else
        {
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

    ReactionRuleDescriptorPyfunc(
        stepladder_type_rrdescriptor stepladder, pyfunc_type pyfunc, boost::shared_ptr<PyObjectHandler> py_handler)
        : stepladder_(stepladder), pyfunc_(pyfunc), pyobject_handler_(py_handler)
    {
        this->pyobject_handler_->inc_ref(this->pyfunc_);
    }

    virtual ~ReactionRuleDescriptorPyfunc()
    {
        this->pyobject_handler_->dec_ref(this->pyfunc_);
    }

    bool is_available() const
    {
        return (pyfunc_ != 0 && stepladder_ != 0);
    }

    virtual Real propensity(const std::vector<Real> &r, const std::vector<Real> &p, Real time) const
    {
        if (stepladder_ == NULL)
        {
            throw IllegalState("stepladder is not registered");
        }
        if (pyfunc_ == NULL)
        {
            throw IllegalState("pyfunc is not registered");
        }
        Real ret = stepladder_(this->pyfunc_, r, p, time);
        return ret;
    }

private:

    pyfunc_type pyfunc_;
    stepladder_type_rrdescriptor stepladder_;
    boost::shared_ptr<PyObjectHandler> pyobject_handler_;
};

} // ecell4

#endif /* ECELL4_REACTION_RULE_DESCRIPTOR_HPP */

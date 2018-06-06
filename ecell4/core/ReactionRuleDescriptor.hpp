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

#include <Python.h>


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

    ReactionRuleDescriptor(const std::string& name)
        : name_(name)
    {
        ;
    }

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

    void set_name(const std::string& name)
    {
        name_ = name;
    }

    const std::string& as_string() const
    {
        return name_;
    }

    // Accessor of coefficients;
    const reaction_coefficient_list_type &reactant_coefficients(void) const
    {
        return this->reactant_coefficients_;
    }

    const reaction_coefficient_list_type &product_coefficients(void) const
    {
        return this->product_coefficients_;
    }

    void set_reactant_coefficient(const std::size_t num, const Real new_coeff)
    {
        this->reactant_coefficients_[num] = new_coeff;
    }

    void set_product_coefficient(const std::size_t num, const Real new_coeff)
    {
        this->product_coefficients_[num] = new_coeff;
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
    std::string name_;
};

class ReactionRuleDescriptorMassAction
    : public ReactionRuleDescriptor
{
public:

    typedef ReactionRuleDescriptor base_type;

public:

    ReactionRuleDescriptorMassAction(const Real k)
        : base_type(""), k_(k)
    {
        ;
    }


    const Real k() const
    {
        return k_;
    }

    void set_k(const Real k)
    {
        k_ = k;
    }

    virtual Real propensity(const std::vector<Real>& reactants, const std::vector<Real>& products, Real t) const
    {
        Real ret = k_;
        for (std::vector<Real>::const_iterator i(reactants.begin()); i != reactants.end(); ++i)
        {
            ret *= (*i);
        }
        return ret;
    }

private:

    Real k_;
};

class ReactionRuleDescriptorCPPfunc
    : public ReactionRuleDescriptor
{
public:

    typedef ReactionRuleDescriptor base_type;
    typedef Real (*func_type)(const std::vector<Real> &r, const std::vector<Real> &p, Real t);

public:

    ReactionRuleDescriptorCPPfunc(func_type pf)
        : base_type(""), pf_(pf)
    {
        ;
    }

    bool is_available() const
    {
        return (this->pf_ != 0);
    }

    func_type get() const
    {
        return this->pf_;
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
    typedef PyObject* pyfunc_type;
    typedef Real (*stepladder_func_type)(pyfunc_type, std::vector<Real>, std::vector<Real>, Real t);

public:

    ReactionRuleDescriptorPyfunc(stepladder_func_type stepladder, pyfunc_type pyfunc, const std::string& name)
        : base_type(name), stepladder_(stepladder), pyfunc_(pyfunc)
    {
        Py_INCREF(this->pyfunc_);
    }

    virtual ~ReactionRuleDescriptorPyfunc()
    {
        Py_DECREF(this->pyfunc_);
    }

    pyfunc_type get() const
    {
        return pyfunc_;
    }

    //XXX: The following implementation doesn't work.
    // bool is_available() const
    // {
    //     return (pyfunc_ != 0 && stepladder_ != 0);
    // }

    virtual Real propensity(const std::vector<Real>& reactants, const std::vector<Real>& products, Real t) const
    {
        // if (stepladder_ == NULL)
        // {
        //     throw IllegalState("stepladder is not registered");
        // }
        // if (pyfunc_ == NULL)
        // {
        //     throw IllegalState("pyfunc is not registered");
        // }
        return this->stepladder_(this->pyfunc_, reactants, products, t);
    }

private:

    pyfunc_type pyfunc_;
    stepladder_func_type stepladder_;
};

} // ecell4

#endif /* ECELL4_REACTION_RULE_DESCRIPTOR_HPP */

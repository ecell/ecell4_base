#ifndef ECELL4_REACTION_RULE_DESCRIPTOR_PYFUNC_HPP
#define ECELL4_REACTION_RULE_DESCRIPTOR_PYFUNC_HPP

#include "Python.h"

#include <stdexcept>
#include <cmath>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#include "types.hpp"
#include "Species.hpp"

#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_io.hpp"

#include "ReactionRuleDescriptor.hpp"


namespace ecell4
{

class ReactionRuleDescriptorPyfunc
    : public ReactionRuleDescriptor
{
public:

    typedef ReactionRuleDescriptor base_type;
    typedef base_type::state_container_type state_container_type;
    typedef base_type::coefficient_container_type coefficient_container_type;
    typedef PyObject* pyfunc_type;
    typedef Real (*stepladder_func_type)(pyfunc_type, const state_container_type&, const state_container_type&, Real volume, Real t, const coefficient_container_type&, const coefficient_container_type&);

public:

    ReactionRuleDescriptorPyfunc(stepladder_func_type stepladder, pyfunc_type pyfunc, const std::string& name)
        : base_type(), stepladder_(stepladder), pyfunc_(pyfunc), name_(name)
    {
        Py_INCREF(this->pyfunc_);
    }

    ReactionRuleDescriptorPyfunc(stepladder_func_type stepladder, pyfunc_type pyfunc, const std::string& name, const coefficient_container_type &reactant_coefficients, const coefficient_container_type &product_coefficients)
        : base_type(reactant_coefficients, product_coefficients), stepladder_(stepladder), pyfunc_(pyfunc), name_(name)
    {
        Py_INCREF(this->pyfunc_);
    }

    virtual ~ReactionRuleDescriptorPyfunc()
    {
        Py_DECREF(this->pyfunc_);
    }

    virtual ReactionRuleDescriptor* clone() const
    {
        return new ReactionRuleDescriptorPyfunc(stepladder_, pyfunc_, as_string(), reactant_coefficients(), product_coefficients());
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

    virtual Real propensity(const state_container_type& reactants, const state_container_type& products, Real volume, Real t) const
    {
        // if (stepladder_ == NULL)
        // {
        //     throw IllegalState("stepladder is not registered");
        // }
        // if (pyfunc_ == NULL)
        // {
        //     throw IllegalState("pyfunc is not registered");
        // }
        return this->stepladder_(this->pyfunc_, reactants, products, volume, t, this->reactant_coefficients(), this->product_coefficients());
    }

    void set_name(const std::string& name)
    {
        name_ = name;
    }

    const std::string& as_string() const
    {
        return name_;
    }

private:

    stepladder_func_type stepladder_;
    pyfunc_type pyfunc_;
    std::string name_;
};

} // ecell4

#endif /* ECELL4_REACTION_RULE_DESCRIPTOR_PYFUNC_HPP */

#ifndef ECELL4_PYTHON_API_REACTION_RULE_DESCRIPTOR_HPP
#define ECELL4_PYTHON_API_REACTION_RULE_DESCRIPTOR_HPP

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <ecell4/core/ReactionRuleDescriptor.hpp>

namespace py = pybind11;

namespace ecell4
{

namespace python_api
{

    template<class Base = ecell4::ReactionRuleDescriptor>
    class PyReactionRuleDescriptor: public Base
    {
    public:
        using Base::Base;
        using state_container_type = ReactionRuleDescriptor::state_container_type;

        Real propensity(const state_container_type& reactants, const state_container_type& products, Real volume, Real t) const
        {
            PYBIND11_OVERLOAD(Real, Base, propensity, reactants, products, volume, t);
        }
    };

    class ReactionRuleDescriptorPyfunc
        : public ReactionRuleDescriptor
    {
    public:
        using base_type = ReactionRuleDescriptor;
        using callback_t = py::object;

        ReactionRuleDescriptorPyfunc(const callback_t& callback, const std::string& name)
            : base_type(), callback_(callback), name_(name)
        {
        }

        ReactionRuleDescriptorPyfunc(const callback_t& callback, const std::string& name,
                const coefficient_container_type& reactant_coefficients,
                const coefficient_container_type& product_coefficients)
            : base_type(reactant_coefficients, product_coefficients), callback_(callback), name_(name)
        {
        }

        Real propensity(const state_container_type& reactants, const state_container_type& products, Real volume, Real t) const override
        {
            return callback_(reactants, products, volume, t, reactant_coefficients(), product_coefficients()).cast<Real>();
        }

        callback_t get() const
        {
            return callback_;
        }

        void set_name(const std::string& name)
        {
            name_ = name;
        }

        const std::string& as_string() const
        {
            return name_;
        }

        ReactionRuleDescriptor* clone() const
        {
            return new ReactionRuleDescriptorPyfunc(callback_, name_,
                    reactant_coefficients(), product_coefficients());
        }

    private:
        callback_t callback_;
        std::string name_;
    };

}

}

#endif /* ECELL4_PYTHON_API_REACTION_RULE_DESCRIPTOR_HPP */

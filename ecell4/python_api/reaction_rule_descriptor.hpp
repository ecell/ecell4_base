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
        using callback_t = std::function<Real (
                const state_container_type&, const state_container_type&,
                Real, Real,
                const coefficient_container_type&, const coefficient_container_type&)>;
        ReactionRuleDescriptorPyfunc(callback_t callback, const std::string& name)
            : base_type(), callback_(callback), name_(name)
        {
        }

        ReactionRuleDescriptorPyfunc(callback_t callback, const std::string& name,
                const coefficient_container_type& reactant_coefficients,
                const coefficient_container_type& product_coefficients)
            : base_type(reactant_coefficients, product_coefficients), callback_(callback), name_(name)
        {
        }

        Real propensity(const state_container_type& reactants, const state_container_type& products, Real volume, Real t) const override
        {
            return callback_(reactants, products, volume, t, reactant_coefficients(), product_coefficients());
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

    private:
        callback_t callback_;
        std::string name_;
    };

    void define_reaction_rule_descriptor(py::module& m)
    {
        py::class_<ReactionRuleDescriptor, PyReactionRuleDescriptor<>,
            boost::shared_ptr<ReactionRuleDescriptor>>(m, "ReactionRuleDescriptor")
            .def("propensity", &ReactionRuleDescriptor::propensity)
            .def("reactant_coefficients", &ReactionRuleDescriptor::reactant_coefficients)
            .def("product_coefficients", &ReactionRuleDescriptor::product_coefficients)
            .def("set_reactant_coefficient", &ReactionRuleDescriptor::set_reactant_coefficient)
            .def("set_product_coefficient", &ReactionRuleDescriptor::set_product_coefficient)
            .def("set_reactant_coefficients", &ReactionRuleDescriptor::set_reactant_coefficients)
            .def("set_product_coefficients", &ReactionRuleDescriptor::set_product_coefficients);

        py::class_<ReactionRuleDescriptorMassAction, ReactionRuleDescriptor,
            PyReactionRuleDescriptor<ReactionRuleDescriptorMassAction>,
            boost::shared_ptr<ReactionRuleDescriptorMassAction>>(m, "ReactionRuleDescriptorMassAction")
            .def(py::init<const Real>())
            .def(py::init<const Quantity<Real>&>())
            .def("k", &ReactionRuleDescriptorMassAction::k)
            .def("get_k", &ReactionRuleDescriptorMassAction::get_k)
            .def("set_k", (void (ReactionRuleDescriptorMassAction::*)(const Real)) &ReactionRuleDescriptorMassAction::set_k)
            .def("set_k", (void (ReactionRuleDescriptorMassAction::*)(const Quantity<Real>&)) &ReactionRuleDescriptorMassAction::set_k)
            .def(py::pickle(
                [](const ReactionRuleDescriptorMassAction& self)
                {
                    return py::make_tuple(self.reactant_coefficients(), self.product_coefficients(), self.get_k());
                },
                [](py::tuple t)
                {
                    if (t.size() != 3)
                        throw std::runtime_error("Invalid state");
                    ReactionRuleDescriptorMassAction ret(t[2].cast<Quantity<Real>>());
                    ret.set_reactant_coefficients(t[0].cast<ReactionRuleDescriptor::coefficient_container_type>());
                    ret.set_product_coefficients(t[1].cast<ReactionRuleDescriptor::coefficient_container_type>());
                    return ret;
                }
            ));

        py::class_<ReactionRuleDescriptorPyfunc, ReactionRuleDescriptor,
            PyReactionRuleDescriptor<ReactionRuleDescriptorPyfunc>,
            boost::shared_ptr<ReactionRuleDescriptorPyfunc>>(m, "ReactionRuleDescriptorPyfunc")
            .def(py::init<ReactionRuleDescriptorPyfunc::callback_t, const std::string&>())
            .def("get", &ReactionRuleDescriptorPyfunc::get)
            .def("set_name", &ReactionRuleDescriptorPyfunc::set_name)
            .def("as_string", &ReactionRuleDescriptorPyfunc::as_string)
            .def(py::pickle(
                [](const ReactionRuleDescriptorPyfunc& self)
                {
                    return py::make_tuple(self.get(), self.as_string(), self.reactant_coefficients(), self.product_coefficients());
                },
                [](py::tuple t)
                {
                    if (t.size() != 4)
                        throw std::runtime_error("Invalid state");
                    return ReactionRuleDescriptorPyfunc(
                        t[0].cast<ReactionRuleDescriptorPyfunc::callback_t>(),
                        t[1].cast<std::string>(),
                        t[2].cast<ReactionRuleDescriptorPyfunc::coefficient_container_type>(),
                        t[3].cast<ReactionRuleDescriptorPyfunc::coefficient_container_type>()
                    );
                }
            ));
    }
}

}

#endif /* ECELL4_PYTHON_API_REACTION_RULE_DESCRIPTOR_HPP */

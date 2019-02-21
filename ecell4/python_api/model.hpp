#ifndef ECELL4_PYTHON_API_MODEL_HPP
#define ECELL4_PYTHON_API_MODEL_HPP

#include <pybind11/pybind11.h>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/NetfreeModel.hpp>

namespace py = pybind11;

namespace ecell4
{

namespace python_api
{

    template<class Base = ecell4::Model>
    class PyModel: public Base
    {
    public:
        using Base::Base;

        std::vector<ReactionRule> query_reaction_rules(const Species& sp) const override
        {
            PYBIND11_OVERLOAD_PURE(std::vector<ReactionRule>, Base, query_reaction_rules, sp);
        }

        std::vector<ReactionRule> query_reaction_rules(const Species& sp1, const Species& sp2) const override
        {
            PYBIND11_OVERLOAD_PURE(std::vector<ReactionRule>, Base, query_reaction_rules, sp1, sp2);
        }

        bool update_species_attribute(const Species& sp) override
        {
            PYBIND11_OVERLOAD(bool, Base, update_species_attribute, sp);
        }

        void add_species_attribute(const Species& sp) override
        {
            PYBIND11_OVERLOAD(void, Base, add_species_attribute, sp);
        }

        bool has_species_attribute(const Species& sp) const override
        {
            PYBIND11_OVERLOAD(bool, Base, has_species_attribute, sp);
        }

        void remove_species_attribute(const Species& sp) override
        {
            PYBIND11_OVERLOAD(void, Base, remove_species_attribute, sp);
        }

        Species apply_species_attributes(const Species& sp) const override
        {
            PYBIND11_OVERLOAD(Species, Base, apply_species_attributes, sp);
        }

        void add_reaction_rule(const ReactionRule& rr) override
        {
            PYBIND11_OVERLOAD(void, Base, add_reaction_rule, rr);
        }

        void remove_reaction_rule(const ReactionRule& rr) override
        {
            PYBIND11_OVERLOAD(void, Base, remove_reaction_rule, rr);
        }

        bool has_reaction_rule(const ReactionRule& rr) const override
        {
            PYBIND11_OVERLOAD(bool, Base, has_reaction_rule, rr);
        }

        const Model::reaction_rule_container_type& reaction_rules() const override
        {
            PYBIND11_OVERLOAD_PURE(const Model::reaction_rule_container_type&, Base, reaction_rules,);
        }

        const Model::species_container_type& species_attributes() const override
        {
            PYBIND11_OVERLOAD_PURE(const Model::species_container_type&, Base, species_attributes,);
        }

        boost::shared_ptr<Model> expand(
            const std::vector<Species>& sp, const Integer max_itr,
            const std::map<Species, Integer>& max_stoich) const override
        {
            PYBIND11_OVERLOAD_PURE(boost::shared_ptr<Model>, Base, expand, sp, max_itr, max_stoich);
        }

        boost::shared_ptr<Model> expand(
            const std::vector<Species>& sp, const Integer max_itr) const override
        {
            PYBIND11_OVERLOAD_PURE(boost::shared_ptr<Model>, Base, expand, sp, max_itr);
        }

        boost::shared_ptr<Model> expand(const std::vector<Species>& sp) const override
        {
            PYBIND11_OVERLOAD_PURE(boost::shared_ptr<Model>, Base, expand, sp);
        }
    };

    template<class Base>
    class PyModelImpl: public PyModel<Base>
    {
    public:
        using PyModel<Base>::PyModel;

        std::vector<ReactionRule> query_reaction_rules(const Species& sp) const override
        {
            PYBIND11_OVERLOAD(std::vector<ReactionRule>, Base, query_reaction_rules, sp);
        }

        std::vector<ReactionRule> query_reaction_rules(const Species& sp1, const Species& sp2) const override
        {
            PYBIND11_OVERLOAD(std::vector<ReactionRule>, Base, query_reaction_rules, sp1, sp2);
        }

        const Model::reaction_rule_container_type& reaction_rules() const override
        {
            PYBIND11_OVERLOAD(const Model::reaction_rule_container_type&, Base, reaction_rules,);
        }

        const Model::species_container_type& species_attributes() const override
        {
            PYBIND11_OVERLOAD(const Model::species_container_type&, Base, species_attributes,);
        }

        boost::shared_ptr<Model> expand(
            const std::vector<Species>& sp, const Integer max_itr,
            const std::map<Species, Integer>& max_stoich) const override
        {
            PYBIND11_OVERLOAD(boost::shared_ptr<Model>, Base, expand, sp, max_itr, max_stoich);
        }

        boost::shared_ptr<Model> expand(
            const std::vector<Species>& sp, const Integer max_itr) const override
        {
            PYBIND11_OVERLOAD(boost::shared_ptr<Model>, Base, expand, sp, max_itr);
        }

        boost::shared_ptr<Model> expand(const std::vector<Species>& sp) const override
        {
            PYBIND11_OVERLOAD(boost::shared_ptr<Model>, Base, expand, sp);
        }
    };

}

}

#endif /* ECELL4_PYTHON_API_MODEL_HPP */

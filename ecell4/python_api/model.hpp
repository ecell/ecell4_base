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

    void define_model(py::module& m)
    {
        py::class_<Model, PyModel<>, boost::shared_ptr<Model>>(m, "Model")
            .def("query_reaction_rules", (std::vector<ReactionRule> (Model::*)(const Species&) const)
                    &Model::query_reaction_rules)
            .def("query_reaction_rules", (std::vector<ReactionRule> (Model::*)(const Species&, const Species&) const)
                    &Model::query_reaction_rules)
            .def("update_species_attribute", &Model::update_species_attribute)
            .def("add_species_attribute", &Model::add_species_attribute)
            .def("has_species_attribute", &Model::has_species_attribute)
            .def("remove_species_attribute", &Model::remove_species_attribute)
            .def("apply_species_attributes", &Model::apply_species_attributes)
            .def("add_reaction_rule", &Model::add_reaction_rule)
            .def("remove_reaction_rule", &Model::remove_reaction_rule)
            .def("has_reaction_rule", &Model::has_reaction_rule)
            .def("reaction_rules", &Model::reaction_rules)
            .def("species_attributes", &Model::species_attributes)
            .def("num_reaction_rules", &Model::num_reaction_rules)
            .def("expand", (boost::shared_ptr<Model> (Model::*)(
                    const std::vector<Species>&, const Integer, const std::map<Species, Integer>&) const) &Model::expand)
            .def("expand", (boost::shared_ptr<Model> (Model::*)(const std::vector<Species>&, const Integer) const) &Model::expand)
            .def("expand", (boost::shared_ptr<Model> (Model::*)(const std::vector<Species>&) const) &Model::expand)
            .def("list_species", &Model::list_species)
            .def("add_species_attributes", &Model::add_species_attributes)
            .def("add_reaction_rules", &Model::add_reaction_rules);

        py::class_<NetworkModel, PyModel<NetworkModel>,
            boost::shared_ptr<NetworkModel>>(m, "NetworkModel")
            .def(py::init<>())
            .def(py::pickle(
                [](const NetworkModel& self)
                {
                    return py::make_tuple(self.species_attributes(), self.reaction_rules());
                },
                [](py::tuple t)
                {
                    if (t.size() != 2)
                        throw std::runtime_error("Invalid state");
                    NetworkModel model;
                    model.add_species_attributes(t[0].cast<Model::species_container_type>());
                    model.add_reaction_rules(t[1].cast<Model::reaction_rule_container_type>());
                    return model;
                }
            ));

        py::class_<NetfreeModel, PyModel<NetfreeModel>,
            boost::shared_ptr<NetfreeModel>>(m, "NetfreeModel")
            .def(py::init<>())
            .def("set_effective", &NetfreeModel::set_effective)
            .def("effective", &NetfreeModel::effective)
            .def(py::pickle(
                [](const NetfreeModel& self)
                {
                    return py::make_tuple(self.species_attributes(), self.reaction_rules());
                },
                [](py::tuple t)
                {
                    if (t.size() != 2)
                        throw std::runtime_error("Invalid state");
                    NetfreeModel model;
                    model.add_species_attributes(t[0].cast<Model::species_container_type>());
                    model.add_reaction_rules(t[1].cast<Model::reaction_rule_container_type>());
                    return model;
                }
            ));
    }
}

}

#endif /* ECELL4_PYTHON_API_MODEL_HPP */

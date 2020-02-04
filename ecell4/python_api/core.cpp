#include "python_api.hpp"

#include <ecell4/core/BDMLWriter.hpp>
#include <ecell4/core/Context.hpp>
#include <ecell4/core/extras.hpp>
#include <ecell4/core/functions.hpp>
#include <ecell4/core/Integer3.hpp>
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/Barycentric.hpp>
#include <ecell4/core/types.hpp>

#include "model.hpp"
#include "observers.hpp"
#include "random_number_generator.hpp"
#include "reaction_rule_descriptor.hpp"
#include "shape.hpp"
#include "world_interface.hpp"
#include "simulator.hpp"

namespace py = pybind11;

namespace ecell4
{

namespace python_api
{

static inline
void define_real3(py::module& m)
{
    py::class_<Real3>(m, "Real3")
        .def(py::init<Real3::value_type, Real3::value_type, Real3::value_type>(),
                py::arg("x"), py::arg("y"), py::arg("z"))
        .def(py::self += py::self)
        .def(py::self + py::self)
        .def(py::self -= py::self)
        .def(py::self - py::self)
        .def(py::self *= Real3::value_type())
        .def(py::self * Real3::value_type())
        .def("__rmul__", [](const Real3& x, Real3::value_type y) { return x * y; }, py::is_operator())
        .def(py::self /= Real3::value_type())
        .def(py::self / Real3::value_type())
        .def("__setitem__",
            [](Real3 &x, std::size_t i, Real3::value_type value)
            {
                if (i >= 3) throw std::out_of_range("");
                x.at(i) = value;
            },
            py::is_operator())
        .def("__getitem__",
            [](const Real3 &x, std::size_t i)
            {
                if (i >= 3) throw std::out_of_range("");
                return x.at(i);
            },
            py::is_operator())
        .def("__abs__", [](const Real3& x) { return abs(x); }, py::is_operator())
        .def("__eq__", [](const Real3& x, const Real3& y)
            {
                return x[0] == y[0] && x[1] == y[1] && x[2] == y[2];
            },
            py::is_operator())
        .def(py::pickle(
            [](const Real3& x)
            {
                return py::make_tuple(x[0], x[1], x[2]);
            },
            [](py::tuple t)
            {
                if (t.size() != 3)
                    throw std::runtime_error("Invalid state");
                return Real3(
                    t[0].cast<Real3::value_type>(),
                    t[1].cast<Real3::value_type>(),
                    t[2].cast<Real3::value_type>());
            }
        ));

    m.def("real3_add", (Real3 (*)(const Real3&, const Real3&)) &add);
    m.def("real3_subtract", (Real3 (*)(const Real3&, const Real3&)) &subtract);
    m.def("real3_divide", (Real3 (*)(const Real3&, const Real3::value_type&)) &divide);
    m.def("real3_multiply", (Real3 (*)(const Real3&, const Real3::value_type&)) &multiply);

    m.def("real3_abs", (Real3 (*)(const Real3&)) abs);
    m.def("real3_dot_product", (Real3::value_type (*)(const Real3&, const Real3&)) &dot_product);
    m.def("cross_product", (Real3 (*)(const Real3&, const Real3&)) &cross_product);

    m.def("real3_length_sq", (Real3::value_type (*)(const Real3&)) &length_sq);
    m.def("real3_length", (Real3::value_type (*)(const Real3&)) &length);

    m.def("length_sq", (Real3::value_type (*)(const Real3&)) &length_sq);
    m.def("length", (Real3::value_type (*)(const Real3&)) &length);
    m.def("dot_product", (Real3::value_type (*)(const Real3&, const Real3&)) &dot_product);

    m.def("ones", &ones);
    m.def("unitx", &unitx);
    m.def("unity", &unity);
    m.def("unitz", &unitz);
}

static inline
void define_integer3(py::module& m)
{
    py::class_<Integer3>(m, "Integer3")
        .def(py::init<Integer3::value_type, Integer3::value_type, Integer3::value_type>(),
                py::arg("col"), py::arg("row"), py::arg("layer"))
        .def_readwrite("col", &Integer3::col)
        .def_readwrite("row", &Integer3::row)
        .def_readwrite("layer", &Integer3::layer)
        .def(py::self == py::self)
        .def(py::self += py::self)
        .def(py::self + py::self)
        .def(py::self -= py::self)
        .def(py::self - py::self)
        .def(py::self *= Integer3::value_type())
        .def("__mul__", [](const Integer3& x, Integer3::value_type y) { return multiply(x, y); }, py::is_operator())
        .def("__rmul__", [](const Integer3& x, Integer3::value_type y) { return multiply(x, y); }, py::is_operator())
        .def("__setitem__", [](Integer3& x, Integer3::size_type i, Integer3::value_type value) { x[i] = value; }, py::is_operator())
        .def("__getitem__", [](const Integer3& x, Integer3::size_type i) { return x[i]; }, py::is_operator())
        .def("__abs__", [](const Integer3& x) { return abs(x); }, py::is_operator())
        .def(py::pickle(
            [](const Integer3& x)
            {
                return py::make_tuple(x.col, x.row, x.layer);
            },
            [](py::tuple t)
            {
                if (t.size() != 3)
                    throw std::runtime_error("Invalid state");
                return Integer3(
                    t[0].cast<Integer3::value_type>(),
                    t[1].cast<Integer3::value_type>(),
                    t[2].cast<Integer3::value_type>());
            }
        ));

    m.def("integer3_add", [](const Integer3& x, const Integer3& y) { return x + y; });
    m.def("integer3_subtract", [](const Integer3& x, const Integer3& y) { return x - y; });
    m.def("integer3_multiply", (Integer3 (*)(const Integer3&, const Integer3::value_type&)) &multiply);

    m.def("integer3_length_sq", (Integer3::value_type (*)(const Integer3&)) &length_sq);
    m.def("integer3_length", (Real (*)(const Integer3&)) &length);

    m.def("integer3_dot_product", (Integer3::value_type (*)(const Integer3&, const Integer3&)) &dot_product);
    m.def("integer3_abs", (Integer3 (*)(const Integer3&)) &abs);

    m.def("length_sq", (Integer3::value_type (*)(const Integer3&)) &length_sq);
    m.def("length", (Real (*)(const Integer3&)) &length);
    m.def("dot_product", (Integer3::value_type (*)(const Integer3&, const Integer3&)) &dot_product);
}

template<typename T>
static inline
py::class_<Quantity<T>> define_quantity(py::module& m, const std::string& name)
{
    using Q = Quantity<T>;
    py::class_<Q> quantity(m, name.c_str());
    quantity
        .def(py::init<const T&, const typename Q::units_type&>(),
                py::arg("magnitude"),
                py::arg("units") = "")
        .def_readwrite("magnitude", &Q::magnitude)
        .def_readwrite("units", &Q::units)
        .def(py::self == py::self)
        .def(py::pickle(
            [](const Q& q)
            {
                return py::make_tuple(q.magnitude, q.units);
            },
            [](py::tuple t)
            {
                if (t.size() != 2)
                    throw std::runtime_error("Invalid state");
                return Q(
                    t[0].cast<T>(),
                    t[1].cast<typename Q::units_type>());
            }
        ));
    return quantity;
}

template<typename T>
static inline
void set_attribute_as(Attribute& attr, const std::pair<std::string, py::object>& key_value)
{
    try
    {
        const T value(key_value.second.cast<T>());
        attr.set(key_value.first, value);
    }
    catch (py::cast_error e)
    {
        ; // do nothing
    }
}

static inline
void define_attribute(py::module& m)
{
    py::class_<Attribute>(m, "Attriubte")
        .def(py::pickle(
            [](const Attribute& self)
            {
                return py::make_tuple(self.values());
            },
            [](py::tuple t)
            {
                if (t.size() != 1)
                {
                    throw std::runtime_error("Invalid state");
                }

                Attribute attr;
                for (const auto& key_value : t[0].cast<std::vector<std::pair<std::string, py::object>>>())
                {
                    set_attribute_as<std::string>(attr, key_value);
                    set_attribute_as<Quantity<Real>>(attr, key_value);
                    set_attribute_as<Quantity<Integer>>(attr, key_value);
                    set_attribute_as<bool>(attr, key_value);
                }
                return attr;
            }
        ));
}

static inline
void define_species(py::module& m)
{
    py::class_<UnitSpecies>(m, "UnitSpecies")
        .def(py::init<>())
        .def(py::init<const std::string&>(), py::arg("name"))
        .def("serial", &UnitSpecies::serial)
        .def("name", &UnitSpecies::name)
        .def("add_site", &UnitSpecies::add_site)
        .def("deserialize", &UnitSpecies::deserialize)
        .def(py::self == py::self)
        .def("__hash__",
            [](const UnitSpecies& self)
            {
                return ECELL4_HASH_STRUCT<UnitSpecies>()(self);
            }
        )
        .def(py::pickle(
            [](const UnitSpecies& self)
            {
                return py::make_tuple(self.serial());
            },
            [](py::tuple t)
            {
                if (t.size() != 1)
                    throw std::runtime_error("Invalid state");
                auto usp = UnitSpecies();
                usp.deserialize(t[0].cast<UnitSpecies::serial_type>());
                return usp;
            }
        ));
    py::class_<Species>(m, "Species")
        .def(py::init<>())
        .def(py::init<const Species::serial_type&>(), py::arg("serial"))
        .def(py::init<const Species::serial_type&, const Real&, const Real&, const std::string, const Integer&>(),
                py::arg("serial"), py::arg("radius"), py::arg("D"),
                py::arg("location") = "",
                py::arg("dimension") = 0)
        .def(py::init<const Species::serial_type&, const Quantity<Real>&, const Quantity<Real>&, const std::string, const Integer&>(),
                py::arg("serial"), py::arg("radius"), py::arg("D"),
                py::arg("location") = "",
                py::arg("dimension") = 0)
        .def("serial", &Species::serial)
        .def("get_attribute", &Species::get_attribute)
        .def("set_attribute", &Species::set_attribute<std::string>)
        .def("set_attribute", &Species::set_attribute<const char*>)
        .def("set_attribute", &Species::set_attribute<bool>)  //XXX: This must be former than Integer's
        .def("set_attribute", &Species::set_attribute<Real>)
        .def("set_attribute", &Species::set_attribute<Integer>)
        .def("set_attribute", &Species::set_attribute<Quantity<Real>>)
        .def("set_attribute", &Species::set_attribute<Quantity<Integer>>)
        .def("remove_attribute", &Species::remove_attribute)
        .def("has_attribute", &Species::has_attribute)
        .def("list_attributes", &Species::list_attributes)
        .def("add_unit", &Species::add_unit)
        .def("count", &Species::count)
        .def("units", &Species::units)
        .def("D", &Species::D)
        .def("radius", &Species::radius)
        .def("location", &Species::location)
        .def("dimension", &Species::dimension)
        .def(py::self == py::self)
        .def(py::self < py::self)
        .def(py::self > py::self)
        .def("__hash__",
            [](const Species& self)
            {
                return ECELL4_HASH_STRUCT<Species>()(self);
            }
        )
        .def(py::pickle(
            [](const Species& species)
            {
                return py::make_tuple(species.serial(), species.attributes());
            },
            [](py::tuple t)
            {
                if (t.size() != 2)
                    throw std::runtime_error("Invalid state");
                Species species(t[0].cast<Species::serial_type>());
                species.set_attributes(t[1].cast<Attribute>());
                return species;
            }
        ));

    m.def("count_species_matches", &count_species_matches);
    m.def("format_species", &format_species);
}

static inline
void define_particle(py::module& m)
{
    py::class_<ParticleID>(m, "ParticleID")
        .def(py::init<>())
        .def(py::init<const ParticleID::value_type>(), py::arg("value"))
        .def("lot", (const ParticleID::lot_type& (ParticleID::*)() const) &ParticleID::lot)
        .def("serial", (const ParticleID::serial_type& (ParticleID::*)() const) &ParticleID::serial)
        .def(py::pickle(
            [](const ParticleID& pid)
            {
                return py::make_tuple(pid.lot(), pid.serial());
            },
            [](py::tuple t)
            {
                if (t.size() != 2)
                    throw std::runtime_error("Invalid state");
                return ParticleID(std::make_pair(
                    t[0].cast<ParticleID::lot_type>(),
                    t[1].cast<ParticleID::serial_type>()
                ));
            }
        ));
    py::class_<Particle>(m, "Particle")
        .def(py::init<const Species&, const Real3&, const Real&, const Real&>(),
                py::arg("sp"), py::arg("pos"), py::arg("radius"), py::arg("D"))
        .def("position", (const Real3& (Particle::*)() const) &Particle::position)
        .def("radius", (const Real& (Particle::*)() const) &Particle::radius)
        .def("D", (const Real& (Particle::*)() const) &Particle::D)
        .def("species", (const Species& (Particle::*)() const) &Particle::species)
        .def(py::pickle(
            [](const Particle& self)
            {
                return py::make_tuple(self.species(), self.position(), self.radius(), self.D());
            },
            [](py::tuple t)
            {
                if (t.size() != 4)
                    throw std::runtime_error("Invalid state");
                return Particle(
                    t[0].cast<Species>(),
                    t[1].cast<Real3>(),
                    t[2].cast<Real>(),
                    t[3].cast<Real>()
                );
            }
        ));
}

static inline
void define_rng(py::module& m)
{
    py::class_<RandomNumberGenerator, PyRandomNumberGenerator<>,
        boost::shared_ptr<RandomNumberGenerator>>(m, "RandomNumberGenerator")
        .def("uniform", &RandomNumberGenerator::uniform)
        .def("uniform", &RandomNumberGenerator::uniform)
        .def("uniform_int", &RandomNumberGenerator::uniform_int)
        .def("gaussian", &RandomNumberGenerator::gaussian,
            py::arg("sigma"), py::arg("mean") = 0.0)
        .def("binomial", &RandomNumberGenerator::binomial)
        .def("seed", (void (RandomNumberGenerator::*)()) &RandomNumberGenerator::seed)
        .def("seed", (void (RandomNumberGenerator::*)(Integer)) &RandomNumberGenerator::seed)
        .def("save", (void (RandomNumberGenerator::*)(const std::string&) const) &RandomNumberGenerator::save)
        .def("load", (void (RandomNumberGenerator::*)(const std::string&)) &RandomNumberGenerator::load);

    py::class_<GSLRandomNumberGenerator, RandomNumberGenerator,
        PyRandomNumberGeneratorImpl<GSLRandomNumberGenerator>,
        boost::shared_ptr<GSLRandomNumberGenerator>>(m, "GSLRandomNumberGenerator")
        .def(py::init<>())
        .def(py::init<const Integer>(), py::arg("seed"))
        .def(py::init<const std::string&>(), py::arg("filename"));
}

static inline
void define_reaction_rule(py::module& m)
{
    using Reactants = ReactionRule::reactant_container_type;
    using Products = ReactionRule::product_container_type;

    py::class_<ReactionRule> reaction_rule(m, "ReactionRule");
    reaction_rule
        .def(py::init<>())
        .def(py::init<const Reactants&, const Products&>(),
                py::arg("reactants"), py::arg("products"))
        .def(py::init<const Reactants&, const Products&, const Real&>(),
                py::arg("reactants"), py::arg("products"), py::arg("k"))
        .def(py::init<const Reactants&, const Products&, const Quantity<Real>&>(),
                py::arg("reactants"), py::arg("products"), py::arg("k"))
        .def("k", &ReactionRule::k)
        .def("set_k", (void (ReactionRule::*)(const Real&)) &ReactionRule::set_k)
        .def("set_k", (void (ReactionRule::*)(const Quantity<Real>&)) &ReactionRule::set_k)
        .def("get_k", &ReactionRule::get_k)
        .def("reactants", &ReactionRule::reactants)
        .def("products", &ReactionRule::products)
        .def("add_reactant", &ReactionRule::add_reactant)
        .def("add_product", &ReactionRule::add_product)
        .def("as_string", &ReactionRule::as_string)
        .def("policy", &ReactionRule::policy)
        .def("set_policy", &ReactionRule::set_policy)
        .def("count", &ReactionRule::count)
        .def("generate", &ReactionRule::generate)
        .def("set_descriptor", &ReactionRule::set_descriptor)
        .def("get_descriptor", &ReactionRule::get_descriptor)
        .def("has_descriptor", &ReactionRule::has_descriptor)
        .def("reset_descriptor", &ReactionRule::reset_descriptor)
        .def("get_attribute", &ReactionRule::get_attribute)
        .def("set_attribute", &ReactionRule::set_attribute<std::string>)
        .def("set_attribute", &ReactionRule::set_attribute<const char*>)
        .def("set_attribute", &ReactionRule::set_attribute<bool>)  //XXX: This must be former than Integer's
        .def("set_attribute", &ReactionRule::set_attribute<Real>)
        .def("set_attribute", &ReactionRule::set_attribute<Integer>)
        .def("set_attribute", &ReactionRule::set_attribute<Quantity<Real>>)
        .def("set_attribute", &ReactionRule::set_attribute<Quantity<Integer>>)
        .def("remove_attribute", &ReactionRule::remove_attribute)
        .def("has_attribute", &ReactionRule::has_attribute)
        .def("list_attributes", &ReactionRule::list_attributes)
        .def(py::self == py::self)
        .def(py::self < py::self)
        .def(py::pickle(
            [](const ReactionRule& self)
            {
                return py::make_tuple(self.reactants(), self.products(), self.get_k(), self.get_descriptor(), self.policy(), self.attributes());
            },
            [](py::tuple t)
            {
                if (t.size() != 6)
                    throw std::runtime_error("Invalid state");
                ReactionRule rr(
                    t[0].cast<Reactants>(),
                    t[1].cast<Products>(),
                    t[2].cast<Quantity<Real> >()
                );
                rr.set_descriptor(t[3].cast<boost::shared_ptr<ReactionRuleDescriptor>>());
                rr.set_policy(t[4].cast<ReactionRule::policy_type>());
                rr.set_attributes(t[5].cast<Attribute>());
                return rr;
            }
        ));

    py::enum_<ReactionRule::policy_type>(reaction_rule, "ReactionRulePolicy")
        .value("STRICT", ReactionRule::policy_type::POLICY_STRICT)
        .value("IMPLICIT", ReactionRule::policy_type::POLICY_IMPLICIT)
        .value("DESTROY", ReactionRule::policy_type::POLICY_DESTROY)
        .export_values()
        .def(py::self | py::self);

    m.def("create_degradation_reaction_rule", &create_degradation_reaction_rule);
    m.def("create_synthesis_reaction_rule", &create_synthesis_reaction_rule);
    m.def("create_unimolecular_reaction_rule", &create_unimolecular_reaction_rule);
    m.def("create_binding_reaction_rule", &create_binding_reaction_rule);
    m.def("create_unbinding_reaction_rule", &create_unbinding_reaction_rule);
}

static inline
void define_model(py::module& m)
{
    py::class_<Model, PyModel<>, boost::shared_ptr<Model>>(m, "Model")
        .def("query_reaction_rules", (std::vector<ReactionRule> (Model::*)(const Species&) const)
                &Model::query_reaction_rules)
        .def("query_reaction_rules", (std::vector<ReactionRule> (Model::*)(const Species&, const Species&) const)
                &Model::query_reaction_rules)
        .def("update_species_attribute", &Model::update_species_attribute)
        .def("add_species_attribute", &Model::add_species_attribute, py::arg("sp"), py::arg("proceed") = false)
        .def("has_species_attribute", &Model::has_species_attribute)
        .def("remove_species_attribute", &Model::remove_species_attribute)
        .def("apply_species_attributes", &Model::apply_species_attributes)
        .def("add_reaction_rule", &Model::add_reaction_rule)
        .def("remove_reaction_rule", &Model::remove_reaction_rule)
        .def("has_reaction_rule", &Model::has_reaction_rule)
        .def("reaction_rules", &Model::reaction_rules)
        .def("species_attributes", &Model::species_attributes)
        .def("species_attributes_proceed", &Model::species_attributes_proceed)
        .def("num_reaction_rules", &Model::num_reaction_rules)
        .def("expand", (boost::shared_ptr<Model> (Model::*)(
                const std::vector<Species>&, const Integer, const std::map<Species, Integer>&) const) &Model::expand)
        .def("expand", (boost::shared_ptr<Model> (Model::*)(const std::vector<Species>&, const Integer) const) &Model::expand)
        .def("expand", (boost::shared_ptr<Model> (Model::*)(const std::vector<Species>&) const) &Model::expand)
        .def("list_species", &Model::list_species)
        .def("add_species_attributes", (void (Model::*)(const std::vector<Species>&)) &Model::add_species_attributes)
        .def("add_species_attributes", (void (Model::*)(const std::vector<std::pair<Species, bool> >&)) &Model::add_species_attributes)
        .def("add_reaction_rules", &Model::add_reaction_rules);

    py::class_<NetworkModel, Model, PyModelImpl<NetworkModel>,
        boost::shared_ptr<NetworkModel>>(m, "NetworkModel")
        .def(py::init<>())
        .def(py::pickle(
            [](const NetworkModel& self)
            {
                return py::make_tuple(
                        self.species_attributes(),
                        self.species_attributes_proceed(),
                        self.reaction_rules());
            },
            [](py::tuple t)
            {
                if (t.size() != 3)
                    throw std::runtime_error("Invalid state");
                NetworkModel model;
                model.add_species_attributes(
                        t[0].cast<Model::species_container_type>(),
                        t[1].cast<std::vector<bool> >());
                model.add_reaction_rules(t[2].cast<Model::reaction_rule_container_type>());
                return model;
            }
        ));

    py::class_<NetfreeModel, Model, PyModelImpl<NetfreeModel>,
        boost::shared_ptr<NetfreeModel>>(m, "NetfreeModel")
        .def(py::init<>())
        .def("set_effective", &NetfreeModel::set_effective)
        .def("effective", &NetfreeModel::effective)
        .def(py::pickle(
            [](const NetfreeModel& self)
            {
                return py::make_tuple(
                        self.species_attributes(),
                        self.species_attributes_proceed(),
                        self.reaction_rules());
            },
            [](py::tuple t)
            {
                if (t.size() != 3)
                    throw std::runtime_error("Invalid state");
                NetfreeModel model;
                model.add_species_attributes(
                        t[0].cast<Model::species_container_type>(),
                        t[1].cast<std::vector<bool> >());
                model.add_reaction_rules(t[2].cast<Model::reaction_rule_container_type>());
                return model;
            }
        ));
}

static inline
void define_world_interface(py::module& m)
{
    py::class_<WorldInterface, PyWorldInterface<>, boost::shared_ptr<WorldInterface>>(m, "WorldInterface")
        .def("t", &WorldInterface::t)
        .def("set_t", &WorldInterface::set_t)
        .def("save", &WorldInterface::save)
        .def("load", &WorldInterface::load)
        .def("volume", &WorldInterface::volume)
        .def("has_species", &WorldInterface::has_species)
        .def("list_species", &WorldInterface::list_species)
        .def("num_molecules", &WorldInterface::num_molecules)
        .def("num_molecules_exact", &WorldInterface::num_molecules_exact)
        .def("get_value", &WorldInterface::get_value)
        .def("get_value_exact", &WorldInterface::get_value_exact)
        .def("edge_lengths", &WorldInterface::edge_lengths)
        .def("num_particles",
            (Integer (WorldInterface::*)() const) &WorldInterface::num_particles)
        .def("num_particles",
            (Integer (WorldInterface::*)(const Species&) const) &WorldInterface::num_particles)
        .def("num_particles_exact", &WorldInterface::num_particles_exact)
        .def("has_particle", &WorldInterface::has_particle)
        .def("get_particle", &WorldInterface::get_particle)
        .def("list_particles",
            (std::vector<std::pair<ParticleID, Particle>> (WorldInterface::*)() const)
            &WorldInterface::list_particles)
        .def("list_particles",
            (std::vector<std::pair<ParticleID, Particle>> (WorldInterface::*)(const Species&) const)
            &WorldInterface::list_particles)
        .def("list_particles_exact", &WorldInterface::list_particles_exact)
        ;
}

static inline
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
        .def(py::init<const Real>(), py::arg("k"))
        .def(py::init<const Quantity<Real>&>(), py::arg("k"))
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
        .def(py::init<ReactionRuleDescriptorPyfunc::callback_t, const std::string&>(),
                py::arg("pyfunc"), py::arg("name"))
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

static inline
void define_observers(py::module& m)
{
    py::class_<NumberLogger>(m, "NumberLogger")
        .def(py::pickle(
            [](const NumberLogger& obj) {
                return py::make_tuple(obj.targets, obj.data, obj.all_species);
                },
            [](py::tuple state) {
                if (state.size() != 3)
                    throw std::runtime_error("Invalid state!");
                auto obj = NumberLogger();
                obj.data = state[1].cast<std::vector<std::vector<Real> > >();
                obj.targets = state[0].cast<std::vector<Species> >();
                obj.all_species = state[2].cast<bool>();
                return obj;
                }
            ));

    py::class_<Observer, PyObserver<>, boost::shared_ptr<Observer>>(m, "Observer")
        .def("next_time", &Observer::next_time)
        .def("reset", &Observer::reset)
        .def("num_steps", &Observer::num_steps);

    py::class_<FixedIntervalPythonHooker, Observer, PyObserver<FixedIntervalPythonHooker>,
        boost::shared_ptr<FixedIntervalPythonHooker>>(m, "FixedIntervalPythonHooker")
        .def(py::init<const Real&, FixedIntervalPythonHooker::callback_t>(),
                py::arg("dt"), py::arg("pyfunc"));

    py::class_<FixedIntervalNumberObserver, Observer, PyObserver<FixedIntervalNumberObserver>,
        boost::shared_ptr<FixedIntervalNumberObserver>>(m, "FixedIntervalNumberObserver")
        .def(py::init<const Real&>(), py::arg("dt"))
        .def(py::init<const Real&, const std::vector<std::string>&>(),
                py::arg("dt"), py::arg("species"))
        .def("data", &FixedIntervalNumberObserver::data)
        .def("targets", &FixedIntervalNumberObserver::targets)
        .def("save", &FixedIntervalNumberObserver::save);

    py::class_<NumberObserver, Observer, PyObserver<NumberObserver>,
        boost::shared_ptr<NumberObserver>>(m, "NumberObserver")
        .def(py::init<>())
        .def(py::init<const std::vector<std::string>&>(), py::arg("species"))
        .def("data", &NumberObserver::data)
        .def("targets", &NumberObserver::targets)
        .def("save", &NumberObserver::save)
        .def(py::pickle(
            [](const NumberObserver& obj) {
                return py::make_tuple(obj.logger(), obj.num_steps());
                },
            [](py::tuple state) {
                if (state.size() != 2)
                    throw std::runtime_error("Invalid state!");
                auto obj = NumberObserver();
                obj.set_logger(state[0].cast<NumberLogger>());
                obj.set_num_steps(state[1].cast<Integer>());
                return obj;
                }
            ));

    py::class_<TimingNumberObserver, Observer, PyObserver<TimingNumberObserver>,
        boost::shared_ptr<TimingNumberObserver>>(m, "TimingNumberObserver")
        .def(py::init<const std::vector<Real>&>(), py::arg("t"))
        .def(py::init<const std::vector<Real>&, const std::vector<std::string>&>(),
                py::arg("t"), py::arg("species"))
        .def("data", &TimingNumberObserver::data)
        .def("targets", &TimingNumberObserver::targets)
        .def("save", &TimingNumberObserver::save)
        .def(py::pickle(
            [](const TimingNumberObserver& obj) {
                return py::make_tuple(obj.logger(), obj.timings(), obj.num_steps(), obj.count());
                },
            [](py::tuple state) {
                if (state.size() != 4)
                    throw std::runtime_error("Invalid state!");
                auto obj = TimingNumberObserver(
                        state[1].cast<std::vector<Real> >(),
                        state[2].cast<Integer>(),
                        state[3].cast<Integer>());
                obj.set_logger(state[0].cast<NumberLogger>());
                // obj.set_num_steps(state[4].cast<Integer>());
                return obj;
                }
            ));

    py::class_<FixedIntervalHDF5Observer, Observer, PyObserver<FixedIntervalHDF5Observer>,
        boost::shared_ptr<FixedIntervalHDF5Observer>>(m, "FixedIntervalHDF5Observer")
        .def(py::init<const Real&, const std::string&>(),
                py::arg("dt"), py::arg("filename"))
        .def("prefix", &FixedIntervalHDF5Observer::prefix)
        .def("filename", (const std::string (FixedIntervalHDF5Observer::*)() const) &FixedIntervalHDF5Observer::filename)
        .def("filename", (const std::string (FixedIntervalHDF5Observer::*)(const Integer) const) &FixedIntervalHDF5Observer::filename)
        .def(py::pickle(
            [](const FixedIntervalHDF5Observer& obj) {
                return py::make_tuple(obj.dt(), obj.t0(), obj.count(), obj.prefix(), obj.num_steps());
                },
            [](py::tuple state) {
                if (state.size() != 5)
                    throw std::runtime_error("Invalid state!");
                auto obj = FixedIntervalHDF5Observer(
                        state[0].cast<Real>(),
                        state[1].cast<Real>(),
                        state[2].cast<Integer>(),
                        state[3].cast<std::string>());
                obj.set_num_steps(state[4].cast<Integer>());
                return obj;
                }
            ));

    py::class_<FixedIntervalCSVObserver, Observer, PyObserver<FixedIntervalCSVObserver>,
        boost::shared_ptr<FixedIntervalCSVObserver>>(m, "FixedIntervalCSVObserver")
        .def(py::init<const Real&, const std::string&>(),
                py::arg("dt"), py::arg("filename"))
        .def(py::init<const Real&, const std::string&, std::vector<std::string>&>(),
                py::arg("dt"), py::arg("filename"), py::arg("species"))
        .def("log", &FixedIntervalCSVObserver::log)
        .def("filename", &FixedIntervalCSVObserver::filename)
        .def("set_header", &FixedIntervalCSVObserver::set_header)
        .def("set_formatter", &FixedIntervalCSVObserver::set_formatter);

    py::class_<CSVObserver, Observer, PyObserver<CSVObserver>, boost::shared_ptr<CSVObserver>>(m, "CSVObserver")
        .def(py::init<const std::string&>(), py::arg("filename"))
        .def(py::init<const std::string&, std::vector<std::string>&>(),
                py::arg("filename"), py::arg("species"))
        .def("log", &CSVObserver::log)
        .def("filename", &CSVObserver::filename)
        .def("set_header", &CSVObserver::set_header)
        .def("set_formatter", &CSVObserver::set_formatter);

    py::class_<FixedIntervalTrajectoryObserver, Observer, PyObserver<FixedIntervalTrajectoryObserver>,
        boost::shared_ptr<FixedIntervalTrajectoryObserver>>(m, "FixedIntervalTrajectoryObserver")
        .def(py::init<const Real&, const std::vector<ParticleID>&, const bool, const Real>(),
                py::arg("dt"), py::arg("pids"),
                py::arg("resolve_boundary") = FixedIntervalTrajectoryObserver::default_resolve_boundary(),
                py::arg("subdt") = FixedIntervalTrajectoryObserver::default_subdt())
        .def(py::init<const Real&, const bool, const Real>(),
                py::arg("dt"),
                py::arg("resolve_boundary") = FixedIntervalTrajectoryObserver::default_resolve_boundary(),
                py::arg("subdt") = FixedIntervalTrajectoryObserver::default_subdt())
        .def("data", &FixedIntervalTrajectoryObserver::data)
        .def("num_tracers", &FixedIntervalTrajectoryObserver::num_tracers)
        .def("t", &FixedIntervalTrajectoryObserver::t);

    py::class_<TimingTrajectoryObserver, Observer, PyObserver<TimingTrajectoryObserver>,
        boost::shared_ptr<TimingTrajectoryObserver>>(m, "TimingTrajectoryObserver")
        .def(py::init<const std::vector<Real>&, const std::vector<ParticleID>&, const bool, const Real>(),
                py::arg("t"), py::arg("pids"),
                py::arg("resolve_boundary") = TimingTrajectoryObserver::default_resolve_boundary(),
                py::arg("subdt") = TimingTrajectoryObserver::default_subdt())
        .def(py::init<const std::vector<Real>&, const bool, const Real>(),
                py::arg("t"),
                py::arg("resolve_boundary") = TimingTrajectoryObserver::default_resolve_boundary(),
                py::arg("subdt") = TimingTrajectoryObserver::default_subdt())
        .def("data", &TimingTrajectoryObserver::data)
        .def("num_tracers", &TimingTrajectoryObserver::num_tracers)
        .def("t", &TimingTrajectoryObserver::t);

    py::class_<FixedIntervalTrackingObserver, Observer, PyObserver<FixedIntervalTrackingObserver>,
        boost::shared_ptr<FixedIntervalTrackingObserver>>(m, "FixedIntervalTrackingObserver")
        .def(py::init<const Real&, const std::vector<Species>&, const bool&, const Real, const Real>(),
                py::arg("dt"), py::arg("species"),
                py::arg("resolve_boundary") = FixedIntervalTrackingObserver::default_resolve_boundary(),
                py::arg("subdt") = FixedIntervalTrackingObserver::default_subdt(),
                py::arg("threshold") = FixedIntervalTrackingObserver::default_threshold())
        .def("data", &FixedIntervalTrackingObserver::data)
        .def("num_tracers", &FixedIntervalTrackingObserver::num_tracers)
        .def("t", &FixedIntervalTrackingObserver::t);

    py::class_<TimeoutObserver, Observer, PyObserver<TimeoutObserver>, boost::shared_ptr<TimeoutObserver>>(m, "TimeoutObserver")
        .def(py::init<>())
        .def(py::init<const Real>(), py::arg("interval"))
        .def("interval", &TimeoutObserver::interval)
        .def("duration", &TimeoutObserver::duration)
        .def("accumulation", &TimeoutObserver::accumulation);
}

static inline
void define_shape(py::module& m)
{
    py::class_<Shape, PyShape<>, boost::shared_ptr<Shape>>(m, "Shape")
        .def("dimension", [](const Shape& self) { return static_cast<Integer>(self.dimension()); })
        .def("is_inside", &Shape::is_inside);

    py::class_<Surface, Shape, PyShapeImpl<Surface>, boost::shared_ptr<Surface>>(m, "Surface")
        .def("root", &Surface::root)
        .def(py::pickle(
            [](const Surface& self)
            {
                return py::make_tuple(self.root());
            },
            [](py::tuple t)
            {
                if (t.size() != 1)
                    throw std::runtime_error("Invalid state");
                return Surface(t[0].cast<const boost::shared_ptr<Shape>&>());
            }
        ));

    py::class_<Union, Shape, PyShapeImpl<Union>, boost::shared_ptr<Union>>(m, "Union")
        .def(py::init<const boost::shared_ptr<Shape>&, const boost::shared_ptr<Shape>&>(),
                py::arg("a"), py::arg("b"))
        .def("surface", &Union::surface)
        .def("one", &Union::one)
        .def("another", &Union::another)
        .def(py::pickle(
            [](const Union& self)
            {
                return py::make_tuple(self.one(), self.another());
            },
            [](py::tuple t)
            {
                if (t.size() != 2)
                    throw std::runtime_error("Invalid state");
                return Union(
                    t[0].cast<const boost::shared_ptr<Shape>&>(),
                    t[1].cast<const boost::shared_ptr<Shape>&>()
                );
            }
        ));

    py::class_<Complement, Shape, PyShapeImpl<Complement>, boost::shared_ptr<Complement>>(m, "Complement")
        .def(py::init<const boost::shared_ptr<Shape>&, const boost::shared_ptr<Shape>&>(),
                py::arg("a"), py::arg("b"))
        .def("surface", &Complement::surface)
        .def("one", &Complement::one)
        .def("another", &Complement::another)
        .def(py::pickle(
            [](const Complement& self)
            {
                return py::make_tuple(self.one(), self.another());
            },
            [](py::tuple t)
            {
                if (t.size() != 2)
                    throw std::runtime_error("Invalid state");
                return Complement(
                    t[0].cast<const boost::shared_ptr<Shape>&>(),
                    t[1].cast<const boost::shared_ptr<Shape>&>()
                );
            }
        ));

    py::class_<AffineTransformation, Shape, PyShapeImpl<AffineTransformation>,
        boost::shared_ptr<AffineTransformation>>(m, "AffineTransformation")
        .def(py::init<>())
        .def(py::init<const boost::shared_ptr<Shape>&>(), py::arg("root"))
        .def(py::init<const boost::shared_ptr<Shape>&, const Real3&, const Real3&, const Real3&, const Real3&>(),
                py::arg("root"), py::arg("first"), py::arg("second"), py::arg("third"), py::arg("shift"))
        .def("translate", &AffineTransformation::translate)
        .def("rescale", &AffineTransformation::rescale)
        .def("xroll", &AffineTransformation::xroll)
        .def("yroll", &AffineTransformation::yroll)
        .def("zroll", &AffineTransformation::zroll)
        .def("surface", &AffineTransformation::surface)
        .def("root", &AffineTransformation::root)
        .def("first", &AffineTransformation::first)
        .def("second", &AffineTransformation::second)
        .def("third", &AffineTransformation::third)
        .def("shift", &AffineTransformation::shift)
        .def(py::pickle(
            [](const AffineTransformation& self)
            {
                return py::make_tuple(self.root(), self.first(), self.second(), self.third(), self.shift());
            },
            [](py::tuple t)
            {
                if (t.size() != 5)
                    throw std::runtime_error("Invalid state");
                return AffineTransformation(
                    t[0].cast<const boost::shared_ptr<Shape>&>(),
                    t[1].cast<const Real3&>(),
                    t[2].cast<const Real3&>(),
                    t[3].cast<const Real3&>(),
                    t[4].cast<const Real3&>()
                );
            }
        ));

    py::class_<Sphere, Shape, PyShapeImpl<Sphere>, boost::shared_ptr<Sphere>>(m, "Sphere")
        .def(py::init<const Real3&, const Real>(), py::arg("center"), py::arg("radius"))
        .def("distance", &Sphere::distance)
        .def("surface", &Sphere::surface)
        .def("center", &Sphere::center)
        .def("radius", &Sphere::radius)
        .def(py::pickle(
            [](const Sphere& self)
            {
                return py::make_tuple(self.center(), self.radius());
            },
            [](py::tuple t)
            {
                if (t.size() != 2)
                    throw std::runtime_error("Invalid state");
                return Sphere(t[0].cast<const Real3&>(), t[1].cast<const Real>());
            }
        ));

    py::class_<SphericalSurface, Shape, PyShapeImpl<SphericalSurface>,
        boost::shared_ptr<SphericalSurface>>(m, "SphericalSurface")
        .def(py::init<const Real3&, const Real>(), py::arg("center"), py::arg("radius"))
        .def("distance", &SphericalSurface::distance)
        .def("inside", &SphericalSurface::inside)
        .def("center", &SphericalSurface::center)
        .def("radius", &SphericalSurface::radius)
        .def(py::pickle(
            [](const SphericalSurface& self)
            {
                return py::make_tuple(self.center(), self.radius());
            },
            [](py::tuple t)
            {
                if (t.size() != 2)
                    throw std::runtime_error("Invalid state");
                return SphericalSurface(t[0].cast<const Real3&>(), t[1].cast<const Real>());
            }
        ));

    py::class_<Cylinder, Shape, PyShapeImpl<Cylinder>, boost::shared_ptr<Cylinder>>(m, "Cylinder")
        .def(py::init<const Real3&, const Real, const Real3&, const Real>(),
                py::arg("center"), py::arg("radius"), py::arg("axis"), py::arg("half_height"))
        .def("distance", &Cylinder::distance)
        .def("surface", &Cylinder::surface)
        .def("center", &Cylinder::center)
        .def("axis", &Cylinder::axis)
        .def("half_height", &Cylinder::half_height)
        .def(py::pickle(
            [](const Cylinder& self)
            {
                return py::make_tuple(self.center(), self.radius(), self.axis(), self.half_height());
            },
            [](py::tuple t)
            {
                if (t.size() != 4)
                    throw std::runtime_error("Invalid state");
                return Cylinder(
                    t[0].cast<const Real3&>(),
                    t[1].cast<const Real>(),
                    t[2].cast<const Real3&>(),
                    t[3].cast<const Real>()
                );
            }
        ));

    py::class_<CylindricalSurface, Shape, PyShapeImpl<CylindricalSurface>,
        boost::shared_ptr<CylindricalSurface>>(m, "CylindricalSurface")
        .def(py::init<const Real3&, const Real, const Real3&, const Real>(),
                py::arg("center"), py::arg("radius"), py::arg("axis"), py::arg("half_height"))
        .def("distance", &CylindricalSurface::distance)
        .def("inside", &CylindricalSurface::inside)
        .def("center", &CylindricalSurface::center)
        .def("radius", &CylindricalSurface::radius)
        .def("axis", &CylindricalSurface::axis)
        .def("half_height", &CylindricalSurface::half_height)
        .def(py::pickle(
            [](const CylindricalSurface& self)
            {
                return py::make_tuple(self.center(), self.radius(), self.axis(), self.half_height());
            },
            [](py::tuple t)
            {
                if (t.size() != 4)
                    throw std::runtime_error("Invalid state");
                return CylindricalSurface(
                    t[0].cast<const Real3&>(),
                    t[1].cast<const Real>(),
                    t[2].cast<const Real3&>(),
                    t[3].cast<const Real>()
                );
            }
        ));

    py::class_<PlanarSurface, Shape, PyShapeImpl<PlanarSurface>,
        boost::shared_ptr<PlanarSurface>>(m, "PlanarSurface")
        .def(py::init<const Real3&, const Real3&, const Real3&>(),
                py::arg("origin"), py::arg("e0"), py::arg("e1"))
        .def("origin", &PlanarSurface::origin)
        .def("e0", &PlanarSurface::e0)
        .def("e1", &PlanarSurface::e1)
        .def("normal", &PlanarSurface::normal)
        .def(py::pickle(
            [](const PlanarSurface& self)
            {
                return py::make_tuple(self.origin(), self.e0(), self.e1());
            },
            [](py::tuple t)
            {
                if (t.size() != 3)
                    throw std::runtime_error("Invalid state");
                return PlanarSurface(
                    t[0].cast<const Real3&>(),
                    t[1].cast<const Real3&>(),
                    t[2].cast<const Real3&>()
                );
            }
        ));

    py::class_<Rod, Shape, PyShapeImpl<Rod>, boost::shared_ptr<Rod>>(m, "Rod")
        .def(py::init<const Real&, const Real&, const Real3&>(),
                py::arg("length"), py::arg("radius"),
                py::arg("origin") = Real3())
        .def("distance", &Rod::distance)
        .def("origin", &Rod::origin)
        .def("length", &Rod::length)
        .def("radius", &Rod::radius)
        .def("shift", &Rod::shift)
        .def("surface", &Rod::surface)
        .def(py::pickle(
            [](const Rod& self)
            {
                return py::make_tuple(self.length(), self.radius(), self.origin());
            },
            [](py::tuple t)
            {
                if (t.size() != 3)
                    throw std::runtime_error("Invalid state");
                return Rod(
                    t[0].cast<const Real>(),
                    t[1].cast<const Real>(),
                    t[2].cast<const Real3&>()
                );
            }
        ));

    py::class_<RodSurface, Shape, PyShapeImpl<RodSurface>, boost::shared_ptr<RodSurface>>(m, "RodSurface")
        .def(py::init<const Real&, const Real&, const Real3&>(),
                py::arg("length"), py::arg("radius"),
                py::arg("origin") = Real3())
        .def("distance", &RodSurface::distance)
        .def("origin", &RodSurface::origin)
        .def("length", &RodSurface::length)
        .def("radius", &RodSurface::radius)
        .def("shift", &RodSurface::shift)
        .def(py::pickle(
            [](const RodSurface& self)
            {
                return py::make_tuple(self.length(), self.radius(), self.origin());
            },
            [](py::tuple t)
            {
                if (t.size() != 3)
                    throw std::runtime_error("Invalid state");
                return RodSurface(
                    t[0].cast<const Real>(),
                    t[1].cast<const Real>(),
                    t[2].cast<const Real3&>()
                );
            }
        ));

    py::class_<AABB, Shape, PyShapeImpl<AABB>, boost::shared_ptr<AABB>>(m, "AABB")
        .def(py::init<const Real3&, const Real3&>(), py::arg("lower"), py::arg("upper"))
        .def("distance", &AABB::distance)
        .def("upper", (const Real3& (AABB::*)() const) &AABB::upper)
        .def("lower", (const Real3& (AABB::*)() const) &AABB::lower)
        .def("upper", (Real3& (AABB::*)()) &AABB::upper)
        .def("lower", (Real3& (AABB::*)()) &AABB::lower)
        .def("surface", &AABB::surface)
        .def(py::pickle(
            [](const AABB& self)
            {
                return py::make_tuple(self.lower(), self.upper());
            },
            [](py::tuple t)
            {
                if (t.size() != 2)
                    throw std::runtime_error("Invalid state");
                return AABB(
                    t[0].cast<const Real3&>(),
                    t[1].cast<const Real3&>()
                );
            }
        ));

    py::class_<MeshSurface, Shape, PyShapeImpl<MeshSurface>, boost::shared_ptr<MeshSurface>>(m, "MeshSurface")
        .def(py::init<const std::string, const Real3&>(), py::arg("filename"), py::arg("edge_lengths"))
        .def("filename", &MeshSurface::filename)
        .def("edge_lengths", &MeshSurface::edge_lengths)
        .def(py::pickle(
            [](const MeshSurface& self)
            {
                return py::make_tuple(self.filename(), self.edge_lengths());
            },
            [](py::tuple t)
            {
                if (t.size() != 2)
                    throw std::runtime_error("Invalid state");
                return MeshSurface(
                    t[0].cast<const std::string>(),
                    t[1].cast<const Real3&>()
                );
            }
        ));

    m.def("create_x_plane",
        [](Real x)
        {
            return PlanarSurface(Real3(x, 0, 0), Real3(0, 1, 0), Real3(0, 0, 1));
        });

    m.def("create_y_plane",
        [](Real y)
        {
            return PlanarSurface(Real3(0, y, 0), Real3(1, 0, 0), Real3(0, 0, 1));
        });

    m.def("create_z_plane",
        [](Real z)
        {
            return PlanarSurface(Real3(0, 0, z), Real3(1, 0, 0), Real3(0, 1, 0));
        });

    // =======================================================================
    // sgfrd polygon related stuff

    py::class_<Triangle, Shape, PyShapeImpl<Triangle>, boost::shared_ptr<Triangle>>(m, "Triangle")
        .def(py::init<const Real3&, const Real3&, const Real3&>(), py::arg("v1"), py::arg("v2"), py::arg("v3"))
        .def("normal",    &Triangle::normal)
        .def("area",      &Triangle::area)
        .def("vertex_at", &Triangle::vertex_at)
        .def("vertices",  &Triangle::vertices)
        .def(py::pickle(
            [](const Triangle& self)
            {
                return py::make_tuple(self.vertex_at(0), self.vertex_at(1), self.vertex_at(2));
            },
            [](py::tuple t)
            {
                if (t.size() != 2)
                {
                    throw std::runtime_error("Invalid state");
                }
                return Triangle(t[0].cast<const Real3&>(), t[1].cast<const Real3&>(), t[2].cast<const Real3&>());
            }
        ));

    py::class_<Barycentric>(m, "Barycentric")
        .def(py::init<const Real, const Real, const Real>())
        .def("__setitem__", [](Barycentric& x, Barycentric::size_type i, Barycentric::value_type value) { x.at(i) = value; }, py::is_operator())
        .def("__getitem__", [](const Barycentric& x, Barycentric::size_type i) { return x.at(i); }, py::is_operator())
        ;

    py::class_<Polygon::FaceID>(m, "FaceID");

    py::class_<Polygon, Shape, PyShapeImpl<Polygon>, boost::shared_ptr<Polygon>>(m, "Polygon")
        .def(py::init<const Real3&, const Integer3&>(), py::arg("edge_lengths"), py::arg("matrix_sizes"))
        .def(py::init<const Real3&, const std::vector<Triangle>&>(), py::arg("edge_lengths"), py::arg("triangles"))
        .def("reset", &Polygon::reset)
        .def("triangles", &Polygon::triangles);

    py::enum_<ecell4::STLFormat>(m, "STLFormat", py::arithmetic())
        .value("Ascii",  ecell4::STLFormat::Ascii)
        .value("Binary", ecell4::STLFormat::Binary)
        .export_values();

    m.def("read_polygon",  &ecell4::read_polygon);
    m.def("write_polygon", &ecell4::write_polygon);
}



static inline
void define_simulator(py::module& m)
{
    py::class_<Simulator, PySimulator<>, boost::shared_ptr<Simulator>>(m, "Simulator")
        .def("initialize", &Simulator::initialize)
        .def("t", &Simulator::t)
        .def("dt", &Simulator::dt)
        .def("set_dt", &Simulator::set_dt)
        .def("num_steps", &Simulator::num_steps)
        .def("step", (void (Simulator::*)()) &Simulator::step)
        .def("step", (bool (Simulator::*)(const Real&)) &Simulator::step)
        .def("check_reaction", &Simulator::check_reaction)
        .def("next_time", &Simulator::next_time);
}

void setup_module(py::module& m)
{
    define_real3(m);
    define_integer3(m);
    const auto quantity_real = define_quantity<Real>(m, "Quantity_Real");
    define_quantity<Integer>(m, "Quantity_Integer");
    m.attr("Quantity") = quantity_real;
    define_attribute(m);
    define_species(m);
    define_particle(m);
    define_rng(m);
    define_reaction_rule(m);
    define_model(m);
    define_world_interface(m);
    define_reaction_rule_descriptor(m);
    define_observers(m);
    define_shape(m);
    define_simulator(m);

    m.def("load_version_information", (std::string (*)(const std::string&)) &extras::load_version_information);
    m.def("get_dimension_from_model", &extras::get_dimension_from_model);
    m.def("get_stoichiometry", &extras::get_stoichiometry);
    m.def("cbrt", &ecell4::cbrt);
    m.attr("N_A") = N_A;
    m.attr("epsilon") = epsilon;
    m.def("_save_bd5", &save_bd5);
}

}

}


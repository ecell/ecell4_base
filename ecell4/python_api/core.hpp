#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/Context.hpp>

namespace py = pybind11;

namespace {
    using namespace ecell4;

    void define_real3(py::module& m)
    {
        py::class_<Real3>(m, "Real3")
            .def(py::init<Real3::value_type, Real3::value_type, Real3::value_type>())
            .def(py::self += py::self)
            .def(py::self + py::self)
            .def(py::self -= py::self)
            .def(py::self - py::self)
            .def(py::self *= Real3::value_type())
            .def(py::self * Real3::value_type())
            .def("__mul__", [](Real3::value_type y, const Real3& x) { return x * y; })
            .def(py::self /= Real3::value_type())
            .def(py::self / Real3::value_type())
            .def("__setitem__", [](Real3 &x, std::size_t i, Real3::value_type value) { x[i] = value; })
            .def("__getitem__", [](const Real3 &x, std::size_t i) { return x[i]; })
            .def("__abs__", [](const Real3& x) { return abs(x); })
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

    void define_integer3(py::module& m)
    {
        py::class_<Integer3>(m, "Integer3")
            .def(py::init<Integer3::value_type, Integer3::value_type, Integer3::value_type>())
            .def_readwrite("col", &Integer3::col)
            .def_readwrite("row", &Integer3::row)
            .def_readwrite("layer", &Integer3::layer)
            .def(py::self += py::self)
            .def(py::self + py::self)
            .def(py::self -= py::self)
            .def(py::self - py::self)
            .def(py::self *= Integer3::value_type())
            .def("__mul__", [](const Integer3& x, Integer3::value_type y) { return multiply(x, y); })
            .def("__mul__", [](Integer3::value_type y, const Integer3& x) { return multiply(x, y); })
            .def("__setitem__", [](Integer3& x, Integer3::size_type i, Integer3::value_type value) { x[i] = value; })
            .def("__getitem__", [](const Integer3& x, Integer3::size_type i) { return x[i]; })
            .def("__abs__", [](const Integer3& x) { return abs(x); })
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
    void set_attribute_as(Species& species, const std::pair<std::string, Species::attribute_type>& key_value)
    {
        if (const T* value = boost::get<T>(&key_value.second))
        {
            species.set_attribute(key_value.first, value);
        }
    }

    template<typename T>
    void define_quantity(py::module& m, const std::string& name)
    {
        using Q = Quantity<T>;
        py::class_<Q>(m, name.c_str())
            .def(py::init<>())
            .def(py::init<const T&, const typename Q::units_type&>())
            .def_readwrite("magnitude", &Q::magnitude)
            .def_readwrite("units", &Q::units);
    }

    void define_species(py::module& m)
    {
        py::class_<UnitSpecies>(m, "UnitSpecies")
            .def(py::init<>())
            .def(py::init<const std::string&>())
            .def("__hash__",
                [](const UnitSpecies& self)
                {
                    return ECELL4_HASH_STRUCT<UnitSpecies>()(self);
                }
            )
            .def("serial", &UnitSpecies::serial)
            .def("name", &UnitSpecies::name)
            .def("add_site", &UnitSpecies::add_site)
            .def("deserialize", &UnitSpecies::deserialize)
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
            .def(py::init<const Species::serial_type&>())
            .def(py::init<const Species::serial_type&, const Real&, const Real&>())
            .def(py::init<const Species::serial_type&, const Real&, const Real&, const std::string>())
            .def(py::init<const Species::serial_type&, const Real&, const Real&, const std::string, const Integer&>())
            .def(py::init<const Species::serial_type&, const Quantity<Real>&, const Quantity<Real>&>())
            .def(py::init<const Species::serial_type&, const Quantity<Real>&, const Quantity<Real>&, const std::string>())
            .def(py::init<const Species::serial_type&, const Quantity<Real>&, const Quantity<Real>&, const std::string, const Integer&>())
            .def("__hash__",
                [](const Species& self)
                {
                    return ECELL4_HASH_STRUCT<Species>()(self);
                }
            )
            .def("serial", &Species::serial)
            .def("get_attribute", &Species::get_attribute)
            .def("set_attribute", &Species::set_attribute<std::string>)
            .def("set_attribute", &Species::set_attribute<const char*>)
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
            .def(py::pickle(
                [](const Species& species)
                {
                    return py::make_tuple(species.serial(), species.list_attributes());
                },
                [](py::tuple t)
                {
                    if (t.size() != 2)
                        throw std::runtime_error("Invalid state");
                    Species species(t[0].cast<Species::serial_type>());
                    for (const auto& key_value : t[1].cast<std::vector<std::pair<std::string, Species::attribute_type>>>())
                    {
                        set_attribute_as<std::string>(species, key_value);
                        set_attribute_as<Quantity<Real>>(species, key_value);
                        set_attribute_as<Quantity<Integer>>(species, key_value);
                        set_attribute_as<bool>(species, key_value);
                    }
                    return species;
                }
            ));

        m.def("count_species_matches", &count_species_matches);
        m.def("format_species", &format_species);
    }

    void define_particle(py::module& m)
    {
        py::class_<ParticleID>(m, "ParticleID")
            .def(py::init<>())
            .def(py::init<const ParticleID::value_type>())
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
            .def(py::init<const Species&, const Real3&, const Real&, const Real&>())
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

    void define_rng(py::module& m)
    {
        py::class_<GSLRandomNumberGenerator>(m, "GSLRandomNumberGenerator")
            .def(py::init<>())
            .def(py::init<const Integer>())
            .def(py::init<const std::string&>())
            .def("uniform", &GSLRandomNumberGenerator::uniform)
            .def("uniform_int", &GSLRandomNumberGenerator::uniform_int)
            .def("gaussian", &GSLRandomNumberGenerator::gaussian,
                py::arg("sigma"), py::arg("mean") = 0.0)
            .def("binomial", &GSLRandomNumberGenerator::binomial)
            .def("seed", (void (GSLRandomNumberGenerator::*)()) &GSLRandomNumberGenerator::seed)
            .def("seed", (void (GSLRandomNumberGenerator::*)(Integer)) &GSLRandomNumberGenerator::seed)
            .def("save", (void (GSLRandomNumberGenerator::*)(const std::string&) const) &GSLRandomNumberGenerator::save)
            .def("load", (void (GSLRandomNumberGenerator::*)(const std::string&)) &GSLRandomNumberGenerator::load);
    }

    void define_reaction_rule(py::module& m)
    {
        using Reactants = ReactionRule::reactant_container_type;
        using Products = ReactionRule::product_container_type;

        py::class_<ReactionRule> reaction_rule(m, "ReactionRule");
        reaction_rule
            .def(py::init<>())
            .def(py::init<const ReactionRule&>())
            .def(py::init<const Reactants&, const Products&>())
            .def(py::init<const Reactants&, const Products&, const Real&>())
            .def(py::init<const Reactants&, const Products&, const Quantity<Real>&>())
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
            .def(py::pickle(
                [](const ReactionRule& self)
                {
                    return py::make_tuple(self.reactants(), self.products(), self.k(), self.get_descriptor());
                },
                [](py::tuple t)
                {
                    if (t.size() != 4)
                        throw std::runtime_error("Invalid state");
                    ReactionRule rr(
                        t[0].cast<Reactants>(),
                        t[1].cast<Products>(),
                        t[2].cast<Real>()
                    );
                    rr.set_descriptor(t[3].cast<boost::shared_ptr<ReactionRuleDescriptor>>());
                    return rr;
                }
            ));

        py::enum_<ReactionRule::policy_type>(reaction_rule, "Policy")
            .value("STRICT", ReactionRule::policy_type::STRICT)
            .value("IMPLICIT", ReactionRule::policy_type::IMPLICIT)
            .value("DESTROY", ReactionRule::policy_type::DESTROY)
            .export_values();

        m.def("create_degradation_reaction_rule", &create_degradation_reaction_rule);
        m.def("create_synthesis_reaction_rule", &create_synthesis_reaction_rule);
        m.def("create_unimolecular_reaction_rule", &create_unimolecular_reaction_rule);
        m.def("create_binding_reaction_rule", &create_binding_reaction_rule);
        m.def("create_unbinding_reaction_rule", &create_unbinding_reaction_rule);
    }

    void define_particle_voxel(py::module& m)
    {
        py::class_<ParticleVoxel>(m, "ParticleVoxel")
            .def(py::init<Species, Integer, Real, Real>())
            .def(py::init<Species, Integer, Real, Real, std::string>())
            .def("coordinate", [](const ParticleVoxel &self) { return self.coordinate; })
            .def("D", [](const ParticleVoxel &self) { return self.D; })
            .def("radius", [](const ParticleVoxel &self) { return self.radius; })
            .def("species", [](const ParticleVoxel &self) { return self.species; })
            .def("loc", [](const ParticleVoxel &self) { return self.loc; })
            .def(py::pickle(
                [](const ParticleVoxel &self)
                {
                    return py::make_tuple(self.species, self.coordinate, self.radius, self.D, self.loc);
                },
                [](py::tuple t)
                {
                    if (t.size() != 5)
                        throw std::runtime_error("Invalid state");
                    return ParticleVoxel(
                        t[0].cast<Species>(),
                        t[1].cast<Integer>(),
                        t[2].cast<Real>(),
                        t[3].cast<Real>(),
                        t[4].cast<std::string>()
                    );
                }
            ));
    }

    void setup_module(py::module& m)
    {
        define_real3(m);
        define_integer3(m);
        define_quantity<Real>(m, "Quantity_Real");
        define_quantity<Integer>(m, "Quantity_Integer");
        define_species(m);
        define_particle(m);
        define_rng(m);
        define_reaction_rule(m);
        define_particle_voxel(m);
    }
}

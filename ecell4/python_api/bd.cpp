#include "python_api.hpp"

#include <ecell4/bd/BDFactory.hpp>
#include <ecell4/bd/BDSimulator.hpp>
#include <ecell4/bd/BDWorld.hpp>

#include "simulator.hpp"
#include "simulator_factory.hpp"
#include "world_interface.hpp"

namespace py = pybind11;
using namespace ecell4::bd;

namespace ecell4
{

namespace python_api
{

static inline
void define_bd_factory(py::module& m)
{
    py::class_<BDFactory> factory(m, "BDFactory");
    factory
        .def(py::init<const Integer3&, Real>(),
                py::arg("matrix_sizes") = BDFactory::default_matrix_sizes(),
                py::arg("bd_dt_factor") = BDFactory::default_bd_dt_factor())
        .def("rng", &BDFactory::rng);
    define_factory_functions(factory);

    m.attr("Factory") = factory;
}

static inline
void define_bd_simulator(py::module& m)
{
    py::class_<BDSimulator, Simulator, PySimulator<BDSimulator>,
        std::shared_ptr<BDSimulator>> simulator(m, "BDSimulator");
    simulator
        .def(py::init<std::shared_ptr<BDWorld>, Real>(),
                py::arg("w"), py::arg("bd_dt_factor") = 1e-5)
        .def(py::init<std::shared_ptr<BDWorld>, std::shared_ptr<Model>, Real>(),
                py::arg("w"), py::arg("m"), py::arg("bd_dt_factor") = 1e-5)
        .def("last_reactions", &BDSimulator::last_reactions)
        .def("set_t", &BDSimulator::set_t);
    define_simulator_functions(simulator);

    m.attr("Simulator") = simulator;
}

static inline
void define_bd_world(py::module& m)
{
    py::class_<BDWorld, WorldInterface, PyWorldImpl<BDWorld>,
        std::shared_ptr<BDWorld>> world(m, "BDWorld");
    world
        .def(py::init<const Real3&, const Integer3&>(),
                py::arg("edge_lengths") = Real3(1.0, 1.0, 1.0),
                py::arg("matrix_sizes") = Integer3(1, 1, 1))
        .def(py::init<const Real3&, const Integer3&, std::shared_ptr<RandomNumberGenerator>>(),
                py::arg("edge_lengths"), py::arg("matrix_sizes"), py::arg("rng"))
        .def(py::init<const std::string&>(), py::arg("filename"))

        .def("new_particle",
            (std::pair<std::pair<ParticleID, Particle>, bool> (BDWorld::*)(const Particle&))
            &BDWorld::new_particle)
        .def("new_particle",
            (std::pair<std::pair<ParticleID, Particle>, bool> (BDWorld::*)(const Species&, const Real3&))
            &BDWorld::new_particle)
        .def("update_particle", &BDWorld::update_particle)
        .def("remove_particle", &BDWorld::remove_particle)
        .def("list_particles_within_radius",
            (std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
             (BDWorld::*)(const Real3&, const Real&) const)
            &BDWorld::list_particles_within_radius)
        .def("list_particles_within_radius",
            (std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
             (BDWorld::*)(const Real3&, const Real&, const ParticleID&) const)
            &BDWorld::list_particles_within_radius)
        .def("list_particles_within_radius",
            (std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
             (BDWorld::*)(const Real3&, const Real&, const ParticleID&, const ParticleID&) const)
            &BDWorld::list_particles_within_radius)
        .def("periodic_transpose", &BDWorld::periodic_transpose)
        .def("apply_boundary", &BDWorld::apply_boundary)
        .def("distance_sq", &BDWorld::distance_sq)
        .def("distance", &BDWorld::distance)
        .def("add_molecules",
            (void (BDWorld::*)(const Species&, const Integer&)) &BDWorld::add_molecules)
        .def("add_molecules",
            (void (BDWorld::*)(const Species&, const Integer&, const std::shared_ptr<Shape>)) &BDWorld::add_molecules)
        .def("remove_molecules", &BDWorld::remove_molecules)
        .def("bind_to", &BDWorld::bind_to)
        .def("rng", &BDWorld::rng);

    m.attr("World") = world;
}

static inline
void define_reaction_info(py::module& m)
{
    using container_type = ReactionInfo::container_type;

    py::class_<ReactionInfo>(m, "ReactionInfo")
        .def(py::init<const Real, const container_type, const container_type>(),
                py::arg("t"), py::arg("reactants"), py::arg("products"))
        .def("t", &ReactionInfo::t)
        .def("reactants", &ReactionInfo::reactants)
        .def("products", &ReactionInfo::products)
        .def(py::pickle(
            [](const ReactionInfo& self)
            {
                return py::make_tuple(self.t(), self.reactants(), self.products());
            },
            [](py::tuple t)
            {
                if (t.size() != 3)
                    throw std::runtime_error("Invalid state");
                return ReactionInfo(
                    t[0].cast<Real>(),
                    t[1].cast<container_type>(),
                    t[2].cast<container_type>()
                );
            }
        ));
}

void setup_bd_module(py::module& m)
{
    define_bd_factory(m);
    define_bd_simulator(m);
    define_bd_world(m);
    define_reaction_info(m);
}

}

}

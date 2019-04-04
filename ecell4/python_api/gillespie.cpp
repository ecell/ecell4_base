#include "python_api.hpp"

#include <ecell4/gillespie/GillespieFactory.hpp>
#include <ecell4/gillespie/GillespieSimulator.hpp>
#include <ecell4/gillespie/GillespieWorld.hpp>

#include "simulator.hpp"
#include "simulator_factory.hpp"
#include "world_interface.hpp"

namespace py = pybind11;
using namespace ecell4::gillespie;

namespace ecell4
{

namespace python_api
{

static inline
void define_gillespie_factory(py::module& m)
{
    py::class_<GillespieFactory> factory(m, "GillespieFactory");
    factory
        .def(py::init<>())
        .def("rng", &GillespieFactory::rng);
    define_factory_functions(factory);

    m.attr("Factory") = factory;
}

static inline
void define_gillespie_simulator(py::module& m)
{
    py::class_<GillespieSimulator, Simulator, PySimulator<GillespieSimulator>,
        boost::shared_ptr<GillespieSimulator>> simulator(m, "GillespieSimulator");
    simulator
        .def(py::init<boost::shared_ptr<GillespieWorld>>())
        .def(py::init<boost::shared_ptr<GillespieWorld>, boost::shared_ptr<Model>>())
        .def("last_reactions", &GillespieSimulator::last_reactions)
        .def("set_t", &GillespieSimulator::set_t);
    define_simulator_functions(simulator);
}

static inline
void define_gillespie_world(py::module& m)
{
    py::class_<GillespieWorld, WorldInterface, PyWorldImpl<GillespieWorld>,
        boost::shared_ptr<GillespieWorld>> world(m, "GillespieWorld");
    world
        .def(py::init<>())
        .def(py::init<const Real3&>())
        .def(py::init<const Real3&, boost::shared_ptr<RandomNumberGenerator>>())
        .def(py::init<const std::string>())
        .def("set_value", &GillespieWorld::set_value)
        .def("add_molecules",
            (void (GillespieWorld::*)(const Species&, const Integer&))
            &GillespieWorld::add_molecules)
        .def("add_molecules",
            (void (GillespieWorld::*)(const Species&, const Integer&, const boost::shared_ptr<Shape>))
            &GillespieWorld::add_molecules)
        .def("remove_molecules", &GillespieWorld::remove_molecules)
        .def("new_particle",
            (std::pair<std::pair<ParticleID, Particle>, bool> (GillespieWorld::*)(const Particle&))
            &GillespieWorld::new_particle)
        .def("new_particle",
            (std::pair<std::pair<ParticleID, Particle>, bool> (GillespieWorld::*)(const Species&, const Real3&))
            &GillespieWorld::new_particle)
        .def("bind_to", &GillespieWorld::bind_to)
        .def("rng", &GillespieWorld::rng);

    m.attr("World") = world;
}

static inline
void define_reaction_info(py::module& m)
{
    using container_type = ReactionInfo::container_type;

    py::class_<ReactionInfo>(m, "ReactionInfo")
        .def(py::init<const Real, const container_type, const container_type>())
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

void setup_gillespie_module(py::module& m)
{
    define_gillespie_factory(m);
    define_gillespie_simulator(m);
    define_gillespie_world(m);
    define_reaction_info(m);
}

}

}

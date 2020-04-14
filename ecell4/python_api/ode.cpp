#include "python_api.hpp"

#include <ecell4/ode/ODEFactory.hpp>
#include <ecell4/ode/ODESimulator.hpp>
#include <ecell4/ode/ODEWorld.hpp>

#include "simulator.hpp"
#include "simulator_factory.hpp"
#include "world_interface.hpp"

namespace py = pybind11;
using namespace ecell4::ode;

namespace ecell4
{

namespace python_api
{

static inline
void define_ode_factory(py::module& m)
{
    py::class_<ODEFactory> factory(m, "ODEFactory");
    factory
        .def(py::init<const ODESolverType, const Real, const Real, const Real>(),
            py::arg("solver_type") = ODEFactory::default_solver_type(),
            py::arg("dt") = ODEFactory::default_dt(),
            py::arg("abs_tol") = ODEFactory::default_abs_tol(),
            py::arg("rel_tol") = ODEFactory::default_rel_tol())
        .def("rng", &ODEFactory::rng);
    define_factory_functions(factory);

    m.attr("Factory") = factory;
}

static inline
void define_ode_simulator(py::module& m)
{
    py::class_<ODESimulator, Simulator, PySimulator<ODESimulator>,
        std::shared_ptr<ODESimulator>> simulator(m, "ODESimulator");
    simulator
        .def(py::init<const std::shared_ptr<ODEWorld>&, const ODESolverType>(),
                py::arg("w"),
                py::arg("solver_type") = ODESolverType::ROSENBROCK4_CONTROLLER)
        .def(py::init<const std::shared_ptr<ODEWorld>&, const std::shared_ptr<Model>&, const ODESolverType>(),
                py::arg("w"), py::arg("m"),
                py::arg("solver_type") = ODESolverType::ROSENBROCK4_CONTROLLER)
        .def("set_t", &ODESimulator::set_t)
        .def("absolute_tolerance", &ODESimulator::absolute_tolerance)
        .def("maximum_step_interval", &ODESimulator::maximum_step_interval)
        .def("set_maximum_step_interval", &ODESimulator::set_maximum_step_interval)
        .def("set_relative_tolerance", &ODESimulator::set_relative_tolerance)
        .def("relative_tolerance", &ODESimulator::relative_tolerance)
        .def("set_absolute_tolerance", &ODESimulator::set_absolute_tolerance)
        .def("values", &ODESimulator::values)
        .def("derivatives", &ODESimulator::derivatives)
        .def("jacobian", &ODESimulator::jacobian)
        .def("fluxes", &ODESimulator::fluxes)
        .def("elasticity", &ODESimulator::elasticity)
        .def("stoichiometry", &ODESimulator::stoichiometry);
    define_simulator_functions(simulator);

    m.attr("Simulator") = simulator;
}

static inline
void define_ode_world(py::module& m)
{
    py::class_<ODEWorld, WorldInterface, PyWorldImpl<ODEWorld>,
        std::shared_ptr<ODEWorld>> world(m, "ODEWorld");
    world
        .def(py::init<const Real3&>(), py::arg("edge_lengths") = Real3(1.0, 1.0, 1.0))
        .def(py::init<const std::string&>(), py::arg("filename"))
        .def("set_volume", &ODEWorld::set_volume)
        .def("new_particle",
            (std::pair<std::pair<ParticleID, Particle>, bool>
             (ODEWorld::*)(const Particle&)) &ODEWorld::new_particle)
        .def("new_particle",
            (std::pair<std::pair<ParticleID, Particle>, bool>
             (ODEWorld::*)(const Species&, const Real3&)) &ODEWorld::new_particle)
        .def("add_molecules",
            (void (ODEWorld::*)(const Species&, const Real&)) &ODEWorld::add_molecules)
        .def("add_molecules",
            (void (ODEWorld::*)(const Species&, const Integer&, const std::shared_ptr<Shape>)) &ODEWorld::add_molecules)
        .def("remove_molecules", &ODEWorld::remove_molecules)
        .def("set_value", &ODEWorld::set_value)
        .def("reserve_species", &ODEWorld::reserve_species)
        .def("release_species", &ODEWorld::release_species)
        .def("bind_to", &ODEWorld::bind_to)
        .def("get_values", &ODEWorld::get_values)
        .def("evaluate",
            (std::vector<Real>
             (ODEWorld::*)(const std::vector<ReactionRule>&) const) &ODEWorld::evaluate)
        .def("evaluate",
            (Real
             (ODEWorld::*)(const ReactionRule&) const) &ODEWorld::evaluate);

    m.attr("World") = world;
}

void setup_ode_module(py::module& m)
{
    py::enum_<ODESolverType>(m, "ODESolverType")
        .value("RUNGE_KUTTA_CASH_KARP54", ODESolverType::RUNGE_KUTTA_CASH_KARP54)
        .value("ROSENBROCK4_CONTROLLER", ODESolverType::ROSENBROCK4_CONTROLLER)
        .value("EULER", ODESolverType::EULER)
        .export_values();

    define_ode_factory(m);
    define_ode_simulator(m);
    define_ode_world(m);
}

}

}

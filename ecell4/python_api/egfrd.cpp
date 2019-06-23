#include "python_api.hpp"

#include <ecell4/egfrd/egfrd.hpp>

#include "simulator.hpp"
#include "simulator_factory.hpp"
#include "world_interface.hpp"

namespace py = pybind11;
// using namespace ecell4::egfrd;

namespace ecell4
{

namespace python_api
{

static inline
void define_bd_factory(py::module& m)
{
    using namespace ::ecell4::egfrd;
    using matrix_sizes_type = EGFRDWorld::matrix_sizes_type;

    py::class_<BDFactory> factory(m, "BDFactory");
    factory
        .def(py::init<const matrix_sizes_type, Real, Integer>(),
            py::arg("matrix_sizes") = BDFactory::default_matrix_sizes(),
            py::arg("bd_dt_factor") = BDFactory::default_bd_dt_factor(),
            py::arg("dissociation_retry_moves") = BDFactory::default_dissociation_retry_moves())
        .def("rng", &BDFactory::rng);
    define_factory_functions(factory);
}

static inline
void define_bd_simulator(py::module& m)
{
    using namespace ::ecell4::egfrd;
    using BDSimulator = egfrd::BDSimulator;
    using world_type = BDSimulator::world_type;
    using model_type = BDSimulator::model_type;

    py::class_<BDSimulator, Simulator, PySimulator<BDSimulator>,
        boost::shared_ptr<BDSimulator>> simulator(m, "BDSimulator");
    simulator
        .def(py::init<boost::shared_ptr<world_type>, Real, int>(),
                py::arg("w"),
                py::arg("bd_dt_factor") = 1.0,
                py::arg("dissociation_retry_moves") = 1)
        .def(py::init<boost::shared_ptr<world_type>, boost::shared_ptr<model_type>, Real, int>(),
                py::arg("w"), py::arg("m"),
                py::arg("bd_dt_factor") = 1.0,
                py::arg("dissociation_retry_moves") = 1)
        .def("last_reactions", &BDSimulator::last_reactions)
        .def("set_t", &BDSimulator::set_t)
        .def("dt_factor", &BDSimulator::dt_factor)
        .def("add_potential",
            (void (BDSimulator::*)(const Species&, const Real&)) &BDSimulator::add_potential)
        .def("add_potential",
            (void (BDSimulator::*)(const Species&, const boost::shared_ptr<Shape>&)) &BDSimulator::add_potential)
        .def("add_potential",
            (void (BDSimulator::*)(const Species&, const boost::shared_ptr<Shape>&, const Real&)) &BDSimulator::add_potential);
    define_simulator_functions(simulator);
}

static inline
void define_egfrd_factory(py::module& m)
{
    using namespace ::ecell4::egfrd;
    using matrix_sizes_type = EGFRDWorld::matrix_sizes_type;

    py::class_<EGFRDFactory> factory(m, "EGFRDFactory");
    factory
        .def(py::init<const matrix_sizes_type, Real, Integer, Real>(),
            py::arg("matrix_sizes") = EGFRDFactory::default_matrix_sizes(),
            py::arg("bd_dt_factor") = EGFRDFactory::default_bd_dt_factor(),
            py::arg("dissociation_retry_moves") = EGFRDFactory::default_dissociation_retry_moves(),
            py::arg("user_max_shell_size") = EGFRDFactory::default_user_max_shell_size())
        .def("rng", &EGFRDFactory::rng);
    define_factory_functions(factory);

    m.attr("Factory") = factory;
}

static inline
void define_egfrd_simulator(py::module& m)
{
    using world_type = ::ecell4::egfrd::EGFRDSimulator::world_type;
    using model_type = ::ecell4::egfrd::EGFRDSimulator::model_type;
    using length_type = ::ecell4::egfrd::EGFRDSimulator::length_type;

    py::class_<::ecell4::egfrd::EGFRDSimulator, Simulator,
               PySimulator<::ecell4::egfrd::EGFRDSimulator>,
               boost::shared_ptr<::ecell4::egfrd::EGFRDSimulator>
                   > simulator(m, "EGFRDSimulator");
    simulator
        .def(py::init<boost::shared_ptr<world_type>, Real, int, length_type>(),
                py::arg("w"),
                py::arg("bd_dt_factor") = 1e-5,
                py::arg("dissociation_retry_moves") = 1,
                py::arg("user_max_shell_size") = std::numeric_limits<length_type>::infinity())
        .def(py::init<boost::shared_ptr<world_type>, boost::shared_ptr<model_type>, Real, int, length_type>(),
                py::arg("w"), py::arg("m"),
                py::arg("bd_dt_factor") = 1e-5,
                py::arg("dissociation_retry_moves") = 1,
                py::arg("user_max_shell_size") = std::numeric_limits<length_type>::infinity())
        .def("last_reactions", &::ecell4::egfrd::EGFRDSimulator::last_reactions)
        .def("set_t", &::ecell4::egfrd::EGFRDSimulator::set_t)
        .def("set_paranoiac", &::ecell4::egfrd::EGFRDSimulator::set_paranoiac);
    define_simulator_functions(simulator);

    m.attr("Simulator") = simulator;
}

static inline
void define_egfrd_world(py::module& m)
{
    using namespace ::ecell4::egfrd;
    using position_type = EGFRDWorld::position_type;
    using matrix_sizes_type = EGFRDWorld::matrix_sizes_type;
    using particle_type = EGFRDWorld::particle_type;
    using particle_id_type = EGFRDWorld::particle_id_type;
    using length_type = EGFRDWorld::length_type;
    using rng_type = EGFRDWorld::rng_type;

    py::class_<EGFRDWorld, WorldInterface, PyWorldImpl<EGFRDWorld>,
        boost::shared_ptr<EGFRDWorld>> world(m, "EGFRDWorld");
    world
        .def(py::init<const position_type&, const matrix_sizes_type&>(),
                py::arg("edge_lengths") = position_type(1.0, 1.0, 1.0),
                py::arg("matrix_sizes") = matrix_sizes_type(3, 3, 3))
        .def(py::init<const position_type&, const matrix_sizes_type&, const boost::shared_ptr<rng_type>&>(),
                py::arg("edge_lengths"),
                py::arg("matrix_sizes"),
                py::arg("rng"))
        .def(py::init<const std::string>(), py::arg("filename"))
        .def("new_particle",
            (std::pair<std::pair<particle_id_type, particle_type>, bool> (EGFRDWorld::*)(const particle_type&))
            &EGFRDWorld::new_particle)
        .def("new_particle",
            (std::pair<std::pair<particle_id_type, particle_type>, bool> (EGFRDWorld::*)(const Species&, const position_type&))
            &EGFRDWorld::new_particle)
        .def("matrix_sizes", &EGFRDWorld::matrix_sizes)
        .def("set_value", &EGFRDWorld::set_value)
        .def("update_particle", &EGFRDWorld::update_particle)
        .def("remove_particle", &EGFRDWorld::remove_particle)
        .def("list_particles_within_radius",
            (std::vector<std::pair<std::pair<particle_id_type, particle_type>, length_type>>
             (EGFRDWorld::*)(const position_type&, const length_type&) const)
            &EGFRDWorld::list_particles_within_radius)
        .def("list_particles_within_radius",
            (std::vector<std::pair<std::pair<particle_id_type, particle_type>, length_type>>
             (EGFRDWorld::*)(const position_type&, const length_type&, const particle_id_type&) const)
            &EGFRDWorld::list_particles_within_radius)
        .def("list_particles_within_radius",
            (std::vector<std::pair<std::pair<particle_id_type, particle_type>, length_type>>
             (EGFRDWorld::*)(const position_type&, const length_type&, const particle_id_type&, const particle_id_type&) const)
            &EGFRDWorld::list_particles_within_radius)
        .def("apply_boundary", &EGFRDWorld::apply_boundary)
        .def("distance",
            (length_type (EGFRDWorld::*)(const position_type&, const position_type&) const) &EGFRDWorld::distance)
        .def("add_molecules",
            (void (EGFRDWorld::*)(const Species&, const Integer&)) &EGFRDWorld::add_molecules)
        .def("add_molecules",
            (void (EGFRDWorld::*)(const Species&, const Integer&, const boost::shared_ptr<Shape>)) &EGFRDWorld::add_molecules)
        .def("remove_molecules", &EGFRDWorld::remove_molecules)
        .def("bind_to", &EGFRDWorld::bind_to)
        .def("rng", &EGFRDWorld::rng);

    m.attr("World") = world;
}

static inline
void define_reaction_info(py::module& m)
{
    using namespace ::ecell4::egfrd;
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

void setup_egfrd_module(py::module& m)
{
    define_bd_simulator(m);
    define_bd_factory(m);
    define_egfrd_factory(m);
    define_egfrd_simulator(m);
    define_egfrd_world(m);
    define_reaction_info(m);
}

}

}

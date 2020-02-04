#include "python_api.hpp"

#include <ecell4/sgfrd/SGFRDFactory.hpp>
#include <ecell4/sgfrd/SGFRDSimulator.hpp>
#include <ecell4/sgfrd/SGFRDWorld.hpp>

#include "simulator.hpp"
#include "simulator_factory.hpp"
#include "world_interface.hpp"

namespace py = pybind11;

namespace ecell4
{
namespace python_api
{

// TODO: enable to read polygon from python function

static inline
void define_sgfrd_factory(py::module& m)
{
    using factory_type = ::ecell4::sgfrd::SGFRDFactory;

    py::class_<factory_type> factory(m, "SGFRDFactory");
    factory
        .def(py::init<const Integer3&, Real, Real>(),
             py::arg("matrix_sizes")              = factory_type::default_matrix_sizes(),
             py::arg("bd_dt_factor")              = factory_type::default_bd_dt_factor(),
             py::arg("bd_reaction_length_factor") = factory_type::default_bd_reaction_length_factor())
        .def("rng",     &factory_type::rng)
        .def("polygon", (factory_type& (factory_type::*)(const Real3&, const std::vector<Triangle>&))&factory_type::polygon)
        .def("polygon", (factory_type& (factory_type::*)(const std::string&, const STLFormat))       &factory_type::polygon);
    define_factory_functions(factory);

    m.attr("Factory") = factory;
}

static inline
void define_sgfrd_simulator(py::module& m)
{
    using factory_type   = ::ecell4::sgfrd::SGFRDFactory;
    using world_type     = ::ecell4::sgfrd::SGFRDWorld;
    using simulator_type = ::ecell4::sgfrd::SGFRDSimulator;

    py::class_<simulator_type, Simulator, PySimulator<simulator_type>,
               boost::shared_ptr<simulator_type>
           > simulator(m, "SGFRDSimulator");
    simulator
        .def(py::init<boost::shared_ptr<world_type>, Real, Real>(),
             py::arg("w"),
             py::arg("bd_dt_factor")              = factory_type::default_bd_dt_factor(),
             py::arg("bd_reaction_length_factor") = factory_type::default_bd_reaction_length_factor())
        .def(py::init<boost::shared_ptr<world_type>, boost::shared_ptr<Model>, Real, Real>(),
             py::arg("w"),
             py::arg("m"),
             py::arg("bd_dt_factor")              = factory_type::default_bd_dt_factor(),
             py::arg("bd_reaction_length_factor") = factory_type::default_bd_reaction_length_factor())
        .def("last_reactions", &simulator_type::last_reactions)
        .def("set_t",          &simulator_type::set_t);
    define_simulator_functions(simulator);

    m.attr("Simulator") = simulator;
}

static inline
void define_sgfrd_world(py::module& m)
{
    using world_type = ::ecell4::sgfrd::SGFRDWorld;

    py::class_<world_type, WorldInterface, PyWorldImpl<world_type>,
        boost::shared_ptr<world_type>> world(m, "SGFRDWorld");
    world
        .def(py::init<const Real3&, const Integer3&>(),
                py::arg("edge_lengths") = Real3(1.0, 1.0, 1.0),
                py::arg("matrix_sizes") = Integer3(1, 1, 1))
        .def(py::init<const Real3&, const Integer3&, boost::shared_ptr<RandomNumberGenerator>>(),
                py::arg("edge_lengths"),
                py::arg("matrix_sizes"),
                py::arg("rng"))
        .def(py::init<const Real3&, const Integer3&, const std::string&, const STLFormat>(),
                py::arg("edge_lengths"),
                py::arg("matrix_sizes"),
                py::arg("stl_file"),
                py::arg("stl_format"))
        .def(py::init<const Real3&, const Integer3&, boost::shared_ptr<RandomNumberGenerator>,
                      const std::string&, const STLFormat>(),
                py::arg("edge_lengths"),
                py::arg("matrix_sizes"),
                py::arg("rng"),
                py::arg("stl_file"),
                py::arg("stl_format"))
        .def(py::init<const std::string&>(), py::arg("filename"))
        .def("polygon", &world_type::polygon)
        .def("new_particle",
            (std::pair<std::pair<ParticleID, Particle>, bool> (world_type::*)(const Particle&))
            &world_type::new_particle)
        .def("new_particle",
            (std::pair<std::pair<ParticleID, Particle>, bool> (world_type::*)(const Species&, const Real3&))
            &world_type::new_particle)
        .def("new_particle",
            (std::pair<std::pair<ParticleID, Particle>, bool>
             (world_type::*)(const Species&, const Polygon::FaceID&, const Barycentric&))
            &world_type::new_particle)
        .def("new_particle",
            (std::pair<std::pair<ParticleID, Particle>, bool>
             (world_type::*)(const Species&, const std::pair<Polygon::FaceID, Barycentric>&))
            &world_type::new_particle)
        .def("update_particle",
            (bool (world_type::*)(const ParticleID&, const Particle&))
            &world_type::update_particle)
        .def("remove_particle",
            (void (world_type::*)(const ParticleID&))
            &world_type::remove_particle)
        .def("list_particles_within_radius",
            (std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
             (world_type::*)(const Real3&, const Real&) const)
            &world_type::list_particles_within_radius)
        .def("list_particles_within_radius",
            (std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
             (world_type::*)(const Real3&, const Real&, const ParticleID&) const)
            &world_type::list_particles_within_radius)
        .def("list_particles_within_radius",
            (std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
             (world_type::*)(const Real3&, const Real&, const ParticleID&, const ParticleID&) const)
            &world_type::list_particles_within_radius)
        .def("periodic_transpose", &world_type::periodic_transpose)
        .def("apply_boundary", &world_type::apply_boundary)
        .def("distance_sq",
            (Real (world_type::*)(const Real3&, const Real3&))
            &world_type::distance_sq)
        .def("distance",
            (Real (world_type::*)(const Real3&, const Real3&))
            &world_type::distance)
        .def("add_molecules",
            (void (world_type::*)(const Species&, const Integer&))
            &world_type::add_molecules)
        .def("add_molecules",
            (void (world_type::*)(const Species&, const Integer&, const boost::shared_ptr<Shape>))
            &world_type::add_molecules)
        .def("remove_molecules", &world_type::remove_molecules)

        .def("get_triangle", &world_type::get_triangle)
        .def("get_surface_position", &world_type::get_surface_position)
        .def("list_surface_positions",
            (std::vector<std::pair<ParticleID, std::pair<Polygon::FaceID, Barycentric>>>
             (world_type::*)() const)
            &world_type::list_surface_positions)
        .def("list_surface_positions",
            (std::vector<std::pair<ParticleID, std::pair<Polygon::FaceID, Barycentric>>>
             (world_type::*)(const Species&) const)
            &world_type::list_surface_positions)
        .def("list_surface_positions_exact", &world_type::list_surface_positions_exact)
        .def("bind_to", &world_type::bind_to)
        .def("rng", (boost::shared_ptr<RandomNumberGenerator>& (world_type::*)())
                    &world_type::rng);
    m.attr("World") = world;
}

static inline
void define_reaction_info(py::module& m)
{
    using reaction_info_type = ::ecell4::sgfrd::ReactionInfo;
    using container_type     = reaction_info_type::container_type;

    py::class_<reaction_info_type>(m, "ReactionInfo")
        .def(py::init<const Real, const container_type, const container_type>(),
                py::arg("t"), py::arg("reactants"), py::arg("products"))
        .def("t",         &reaction_info_type::t)
        .def("reactants", &reaction_info_type::reactants)
        .def("products",  &reaction_info_type::products)
        .def(py::pickle(
            [](const reaction_info_type& self)
            {
                return py::make_tuple(self.t(), self.reactants(), self.products());
            },
            [](py::tuple t)
            {
                if (t.size() != 3)
                    throw std::runtime_error("Invalid state");
                return reaction_info_type(
                    t[0].cast<Real>(),
                    t[1].cast<container_type>(),
                    t[2].cast<container_type>()
                );
            }
        ));
}

void setup_sgfrd_module(py::module& m)
{
    define_sgfrd_factory(m);
    define_sgfrd_simulator(m);
    define_sgfrd_world(m);
    define_reaction_info(m);
}

} // python_api
} // ecell4

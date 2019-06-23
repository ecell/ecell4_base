#include "python_api.hpp"

#include <ecell4/meso/MesoscopicFactory.hpp>
#include <ecell4/meso/MesoscopicSimulator.hpp>
#include <ecell4/meso/MesoscopicWorld.hpp>

#include "simulator.hpp"
#include "simulator_factory.hpp"
#include "world_interface.hpp"

namespace py = pybind11;
using namespace ecell4::meso;

namespace ecell4
{

namespace python_api
{

static inline
void define_meso_factory(py::module& m)
{
    py::class_<MesoscopicFactory> factory(m, "MesoscopicFactory");
    factory
        .def(py::init<const Integer3&, const Real>(),
            py::arg("matrix_sizes") = MesoscopicFactory::default_matrix_sizes(),
            py::arg("subvolume_length") = MesoscopicFactory::default_subvolume_length())
        .def("rng", &MesoscopicFactory::rng);
    define_factory_functions(factory);

    m.attr("Factory") = factory;
}

static inline
void define_meso_simulator(py::module& m)
{
    py::class_<MesoscopicSimulator, Simulator, PySimulator<MesoscopicSimulator>,
        boost::shared_ptr<MesoscopicSimulator>> simulator(m, "MesoscopicSimulator");
    simulator
        .def(py::init<boost::shared_ptr<MesoscopicWorld>>(), py::arg("w"))
        .def(py::init<boost::shared_ptr<MesoscopicWorld>, boost::shared_ptr<Model>>(),
                py::arg("w"), py::arg("m"))
        .def("last_reactions", &MesoscopicSimulator::last_reactions)
        .def("set_t", &MesoscopicSimulator::set_t);
    define_simulator_functions(simulator);

    m.attr("Simulator") = simulator;
}

static inline
void define_meso_world(py::module& m)
{
    using coordinate_type = MesoscopicWorld::coordinate_type;

    py::class_<MesoscopicWorld, WorldInterface, PyWorldImpl<MesoscopicWorld>,
        boost::shared_ptr<MesoscopicWorld>> world(m, "MesoscopicWorld");
    world
        .def(py::init<const Real3&>(), py::arg("edge_lengths") = Real3(1.0, 1.0, 1.0))
        .def(py::init<const Real3&, const Integer3&>(),
                py::arg("edge_lengths"), py::arg("matrix_sizes"))
        .def(py::init<const Real3&, const Integer3&, boost::shared_ptr<RandomNumberGenerator>>(),
                py::arg("edge_lengths"), py::arg("matrix_sizes"), py::arg("rng"))
        .def(py::init<const Real3&, const Real>(),
                py::arg("edge_lengths"), py::arg("subvlume_length"))
        .def(py::init<const Real3&, const Real, boost::shared_ptr<RandomNumberGenerator>>(),
                py::arg("edge_lengths"), py::arg("subvlume_length"), py::arg("rng"))
        .def(py::init<const std::string&>(), py::arg("filename"))
        .def("matrix_sizes", &MesoscopicWorld::matrix_sizes)
        .def("subvolume", &MesoscopicWorld::subvolume)
        .def("set_value", &MesoscopicWorld::set_value)
        .def("num_subvolumes",
            (const Integer (MesoscopicWorld::*)() const) &MesoscopicWorld::num_subvolumes)
        .def("num_subvolumes",
            (const Integer (MesoscopicWorld::*)(const Species&) const) &MesoscopicWorld::num_subvolumes)
        .def("subvolume_edge_lengths", &MesoscopicWorld::subvolume_edge_lengths)
        .def("global2coord", &MesoscopicWorld::global2coord)
        .def("coord2global", &MesoscopicWorld::coord2global)
        .def("position2global", &MesoscopicWorld::position2global)
        .def("position2coordinate", &MesoscopicWorld::position2coordinate)
        .def("num_molecules",
            (Integer (MesoscopicWorld::*)(const Species&) const) &MesoscopicWorld::num_molecules)
        .def("num_molecules",
            (Integer (MesoscopicWorld::*)(const Species&, const coordinate_type&) const) &MesoscopicWorld::num_molecules)
        .def("num_molecules",
            (Integer (MesoscopicWorld::*)(const Species&, const Integer3&) const) &MesoscopicWorld::num_molecules)
        .def("add_molecules",
            (void (MesoscopicWorld::*)(const Species&, const Integer&)) &MesoscopicWorld::add_molecules)
        .def("add_molecules",
            (void (MesoscopicWorld::*)(const Species&, const Integer&, const coordinate_type&)) &MesoscopicWorld::add_molecules)
        .def("add_molecules",
            (void (MesoscopicWorld::*)(const Species&, const Integer&, const Integer3&)) &MesoscopicWorld::add_molecules)
        .def("add_molecules",
            (void (MesoscopicWorld::*)(const Species&, const Integer&, const boost::shared_ptr<Shape>)) &MesoscopicWorld::add_molecules)
        .def("remove_molecules",
            (void (MesoscopicWorld::*)(const Species&, const Integer&)) &MesoscopicWorld::remove_molecules)
        .def("remove_molecules",
            (void (MesoscopicWorld::*)(const Species&, const Integer&, const Integer3&)) &MesoscopicWorld::remove_molecules)
        .def("remove_molecules",
            (void (MesoscopicWorld::*)(const Species&, const Integer&, const coordinate_type&)) &MesoscopicWorld::remove_molecules)
        .def("add_structure", &MesoscopicWorld::add_structure)
        .def("get_volume", &MesoscopicWorld::get_volume)
        .def("get_data", &MesoscopicWorld::get_volume)
        .def("has_structure", &MesoscopicWorld::has_structure)
        .def("on_structure",
            (bool (MesoscopicWorld::*)(const Species&, const Integer3&) const)
            &MesoscopicWorld::on_structure)
        .def("check_structure",
            (bool (MesoscopicWorld::*)(const Species&, const Integer3&) const)
            &MesoscopicWorld::check_structure)
        .def("get_occupancy",
            (Real (MesoscopicWorld::*)(const Species&, const coordinate_type&) const)
            &MesoscopicWorld::get_occupancy)
        .def("get_occupancy",
            (Real (MesoscopicWorld::*)(const Species&, const Integer3&) const)
            &MesoscopicWorld::get_occupancy)
        .def("list_coordinates", &MesoscopicWorld::list_coordinates)
        .def("list_coordinates_exact", &MesoscopicWorld::list_coordinates_exact)
        .def("new_particle",
            (std::pair<std::pair<ParticleID, Particle>, bool> (MesoscopicWorld::*)(const Particle&))
            &MesoscopicWorld::new_particle)
        .def("new_particle",
            (std::pair<std::pair<ParticleID, Particle>, bool> (MesoscopicWorld::*)(const Species&, const Real3&))
            &MesoscopicWorld::new_particle)
        .def("bind_to", &MesoscopicWorld::bind_to)
        .def("rng", &MesoscopicWorld::rng);

    m.attr("World") = world;
}

static inline
void define_reaction_info(py::module& m)
{
    using container_type = ReactionInfo::container_type;
    using coordinate_type = ReactionInfo::coordinate_type;

    py::class_<ReactionInfo>(m, "ReactionInfo")
        .def(py::init<const Real, const container_type, const container_type, const coordinate_type>(),
                py::arg("t"), py::arg("reactants"), py::arg("products"), py::arg("coord"))
        .def("t", &ReactionInfo::t)
        .def("coordinate", &ReactionInfo::coordinate)
        .def("reactants", &ReactionInfo::reactants)
        .def("products", &ReactionInfo::products)
        .def(py::pickle(
            [](const ReactionInfo& self)
            {
                return py::make_tuple(self.t(), self.reactants(), self.products(), self.coordinate());
            },
            [](py::tuple t)
            {
                if (t.size() != 4)
                    throw std::runtime_error("Invalid state");
                return ReactionInfo(
                    t[0].cast<Real>(),
                    t[1].cast<container_type>(),
                    t[2].cast<container_type>(),
                    t[3].cast<coordinate_type>()
                );
            }
        ));
}

void setup_meso_module(py::module& m)
{
    define_meso_factory(m);
    define_meso_simulator(m);
    define_meso_world(m);
    define_reaction_info(m);
}

}

}

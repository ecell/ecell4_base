#include "python_api.hpp"

#include <ecell4/spatiocyte/SpatiocyteReactions.hpp>
#include <ecell4/spatiocyte/SpatiocyteFactory.hpp>
#include <ecell4/spatiocyte/SpatiocyteSimulator.hpp>
#include <ecell4/spatiocyte/SpatiocyteWorld.hpp>
#include <ecell4/spatiocyte/Voxel.hpp>

#include "simulator.hpp"
#include "simulator_factory.hpp"
#include "world_interface.hpp"

namespace py = pybind11;
using namespace ecell4::spatiocyte;

namespace ecell4
{

namespace python_api
{

static inline
void define_reaction_info(py::module& m)
{
    using container_type = ReactionInfo::container_type;
    py::class_<ReactionInfo>(m, "ReactionInfo")
        .def(py::init<const Real, const container_type&, const container_type&>())
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

    py::class_<ReactionInfo::Item>(m, "ReactionInfoItem")
        .def_readonly("pid", &ReactionInfo::Item::pid)
        .def_readonly("species", &ReactionInfo::Item::species)
        .def_readonly("voxel", &ReactionInfo::Item::voxel);
}

static inline
void define_spatiocyte_factory(py::module& m)
{
    using world_type = SpatiocyteFactory::world_type;
    using simulator_type = SpatiocyteFactory::simulator_type;

    py::class_<SpatiocyteFactory> factory(m, "SpatiocyteFactory");
    factory
        .def(py::init<>())
        .def(py::init<const Real>())
        .def("rng", &SpatiocyteFactory::rng);
    define_factory_functions(factory);
}

static inline
void define_spatiocyte_simulator(py::module& m)
{
    py::class_<SpatiocyteSimulator, Simulator, PySimulator<SpatiocyteSimulator>>(m, "SpatiocyteSimulator")
        .def(py::init<boost::shared_ptr<SpatiocyteWorld>>())
        .def(py::init<boost::shared_ptr<SpatiocyteWorld>, boost::shared_ptr<Model>>())
        .def("last_reactions", &SpatiocyteSimulator::last_reactions)
        .def("model", &SpatiocyteSimulator::model)
        .def("world", &SpatiocyteSimulator::world)
        .def("run",
            (void (SpatiocyteSimulator::*)(const Real&, const bool)) &SpatiocyteSimulator::run,
            py::arg("duration"), py::arg("is_dirty") = true)
        .def("run",
            (void (SpatiocyteSimulator::*)(const Real&, const boost::shared_ptr<Observer>&, const bool)) &SpatiocyteSimulator::run,
            py::arg("duration"), py::arg("observer"), py::arg("is_dirty") = true)
        .def("run",
            (void (SpatiocyteSimulator::*)(const Real&, std::vector<boost::shared_ptr<Observer>>, const bool)) &SpatiocyteSimulator::run,
            py::arg("duration"), py::arg("observers"), py::arg("is_dirty") = true)
        ;
}

static inline
void define_spatiocyte_world(py::module& m)
{
    py::class_<SpatiocyteWorld, ecell4::WorldInterface, PyWorldImpl<SpatiocyteWorld>,
        boost::shared_ptr<SpatiocyteWorld>>(m, "SpatiocyteWorld")
        .def(py::init<>())
        .def("new_particle", (boost::optional<ParticleID> (SpatiocyteWorld::*)(const Particle&)) &SpatiocyteWorld::new_particle)
        .def("new_particle", (boost::optional<ParticleID> (SpatiocyteWorld::*)(const Species&, const Real3&)) &SpatiocyteWorld::new_particle)
        .def("remove_particle", &SpatiocyteWorld::remove_particle)
        .def("list_structure_particles", &SpatiocyteWorld::list_structure_particles)
        .def("list_non_structure_particles", &SpatiocyteWorld::list_non_structure_particles)
        .def("update_particle", &SpatiocyteWorld::update_particle)
        .def("add_molecules",
                (bool (SpatiocyteWorld::*)(const Species&, const Integer&))
                &SpatiocyteWorld::add_molecules)
        .def("add_molecules",
                (bool (SpatiocyteWorld::*)(const Species&, const Integer&, const boost::shared_ptr<const Shape>))
                &SpatiocyteWorld::add_molecules)
        .def("remove_molecules", &SpatiocyteWorld::remove_molecules)
        .def("voxel_volume", &SpatiocyteWorld::voxel_volume)
        .def("get_volume", &SpatiocyteWorld::get_volume)
        .def("get_voxel", &SpatiocyteWorld::find_voxel)
        .def("get_voxel_at", &SpatiocyteWorld::get_voxel_at)
        .def("set_value", &SpatiocyteWorld::set_value)
        .def("new_voxel", &SpatiocyteWorld::new_voxel)
        .def("new_voxel_structure", &SpatiocyteWorld::new_voxel_structure)
        .def("voxel_radius", &SpatiocyteWorld::voxel_radius)
        .def("size", &SpatiocyteWorld::size)
        .def("shape", &SpatiocyteWorld::shape)
        .def("bind_to", &SpatiocyteWorld::bind_to)
        .def("get_voxel_near_by", &SpatiocyteWorld::position2voxel)
        .def("add_structure", &SpatiocyteWorld::add_structure)
        .def("rng", &SpatiocyteWorld::rng)
        .def_static("calculate_voxel_volume", &SpatiocyteWorld::calculate_voxel_volume)
        .def_static("calculate_hcp_lengths", &SpatiocyteWorld::calculate_hcp_lengths)
        .def_static("calculate_shape", &SpatiocyteWorld::calculate_shape)
        .def_static("calculate_volume", &SpatiocyteWorld::calculate_volume);

    m.def("create_spatiocyte_world_cell_list_impl", &create_spatiocyte_world_cell_list_impl_alias);
    m.def("create_spatiocyte_world_vector_impl", &create_spatiocyte_world_vector_impl_alias);
    m.def("create_spatiocyte_world_square_offlattice_impl", &allocate_spatiocyte_world_square_offlattice_impl);
}

static inline
void define_voxel(py::module& m)
{
    py::class_<Voxel>(m, "Voxel")
        .def("position", &Voxel::position)
        .def("list_neighbors",
            [](const Voxel& self)
            {
                std::vector<Voxel> list;
                for (auto i = 0; i < self.num_neighbors(); ++i)
                {
                    list.push_back(self.get_neighbor(i));
                }
                return list;
            });
}

void setup_spatiocyte_module(py::module& m)
{
    define_reaction_info(m);
    define_spatiocyte_factory(m);
    define_spatiocyte_simulator(m);
    define_spatiocyte_world(m);
    define_voxel(m);
}

}

}


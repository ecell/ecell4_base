from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *


## Cpp_LatticeWorld
#  ecell4::lattice::LatticeWorld
cdef extern from "ecell4/lattice/LatticeWorld.hpp" namespace "ecell4::lattice":
    cdef cppclass Cpp_LatticeWorld "ecell4::lattice::LatticeWorld":
        Cpp_LatticeWorld(
            Cpp_Real3& edge_lengths, const Real& voxel_radius,
            shared_ptr[Cpp_RandomNumberGenerator] rng) except +
        Cpp_LatticeWorld(
            Cpp_Real3& edge_lengths, const Real& voxel_radius) except +
        Cpp_LatticeWorld(Cpp_Real3& edge_lengths) except +
        Cpp_LatticeWorld(string&) except +
        Cpp_LatticeWorld() except +

        void set_t(Real t)
        Real t()
        Cpp_Real3 edge_lengths()
        Real volume()
        Real voxel_volume()
        Cpp_Real3 actual_lengths()
        Real actual_volume()

        pair[pair[Cpp_ParticleID, Cpp_Particle], bool] new_particle(Cpp_Particle& p)
        pair[pair[Cpp_ParticleID, Cpp_Particle], bool] new_particle(Cpp_Species& sp, Cpp_Real3& pos)
        bool remove_particle(Cpp_ParticleID& pid)
        bool remove_voxel(Cpp_ParticleID& pid)
        pair[Cpp_ParticleID, Cpp_Particle] get_particle(Cpp_ParticleID& pid)
        pair[Cpp_ParticleID, Cpp_Voxel] get_voxel(Cpp_ParticleID& pid)

        Integer num_particles()
        Integer num_particles(Cpp_Species& sp)
        Integer num_particles_exact(Cpp_Species& sp)
        Integer num_voxels()
        Integer num_voxels(Cpp_Species& sp)
        Integer num_voxels_exact(Cpp_Species& sp)
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles()
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles(Cpp_Species& sp)
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles_exact(Cpp_Species& sp)
        bool has_particle(Cpp_ParticleID& pid)
        bool update_particle(Cpp_ParticleID& pid, Cpp_Particle& p)
        # vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]] list_particles_within_radius(Cpp_Real3& pos, Real& radius)
        # vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]] list_particles_within_radius(Cpp_Real3& pos, Real& radius, Cpp_ParticleID& ignore)
        # vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]] list_particles_within_radius(Cpp_Real3& pos, Real& radius, Cpp_ParticleID& ignore1, Cpp_ParticleID& ignore2)
        # Cpp_Real3 periodic_transpose(Cpp_Real3& pos1, Cpp_Real3& pos2)
        # Cpp_Real3 apply_boundary(Cpp_Real3& pos)
        # Real distance_sq(Cpp_Real3& pos1, Cpp_Real3& pos2)
        # Real distance(Cpp_Real3& pos1, Cpp_Real3& pos2)
        # # bool has_species(Cpp_Species& sp)
        Integer num_molecules()
        Integer num_molecules(Cpp_Species& sp)
        Integer num_molecules_exact(Cpp_Species& sp)
        void add_molecules(Cpp_Species& sp, Integer num)
        void remove_molecules(Cpp_Species& sp, Integer num)
        # shared_ptr[Cpp_GSLRandomNumberGenerator] rng()
        Integer get_neighbor(Integer, Integer)
        void save(string filename) except +
        void load(string filename)
        pair[pair[Cpp_ParticleID, Cpp_Voxel], bool] new_voxel(Cpp_Voxel& p)
        pair[pair[Cpp_ParticleID, Cpp_Voxel], bool] new_voxel(Cpp_Species& sp, Integer pos)
        vector[pair[Cpp_ParticleID, Cpp_Voxel]] list_voxels()
        vector[pair[Cpp_ParticleID, Cpp_Voxel]] list_voxels(Cpp_Species& sp)
        vector[pair[Cpp_ParticleID, Cpp_Voxel]] list_voxels_exact(Cpp_Species& sp)
        bool update_voxel(Cpp_ParticleID, Cpp_Voxel)
        bool has_voxel(Cpp_ParticleID)
        Real voxel_radius()
        Integer col_size()
        Integer row_size()
        Integer layer_size()
        Integer size()
        Cpp_Integer3 shape()
        void bind_to(shared_ptr[Cpp_Model])
        Cpp_Real3 coordinate2position(Integer)
        Integer position2coordinate(Cpp_Real3)
        shared_ptr[Cpp_RandomNumberGenerator] rng()

        Cpp_Real3 private2position(Integer)
        Integer private2coord(Integer)
        Integer coord2private(Integer)
        Cpp_Integer3 coord2global(Integer)
        Integer global2coord(Cpp_Integer3)
        Cpp_Integer3 private2global(Integer)
        Integer global2private(Cpp_Integer3)
        Cpp_Real3 global2position(Cpp_Integer3)
        Cpp_Integer3 position2global(Cpp_Real3)
        Integer add_structure(Cpp_Species&, shared_ptr[Cpp_Shape])
        void add_molecules(Cpp_Species& sp, Integer num, shared_ptr[Cpp_Shape])

    cdef Cpp_LatticeWorld* create_lattice_world_cell_list_impl_alias(
        Cpp_Real3&, Real, Cpp_Integer3&, shared_ptr[Cpp_RandomNumberGenerator]&)
    cdef Cpp_LatticeWorld* create_lattice_world_vector_impl_alias(
        Cpp_Real3&, Real, shared_ptr[Cpp_RandomNumberGenerator]&)

## LatticeWorld
#  a python wrapper for Cpp_LatticeWorld
cdef class LatticeWorld:
    cdef shared_ptr[Cpp_LatticeWorld]* thisptr

cdef LatticeWorld LatticeWorld_from_Cpp_LatticeWorld(
    shared_ptr[Cpp_LatticeWorld] m)

## Cpp_LatticeSimulator
#  ecell4::lattice::LatticeSimulator
cdef extern from "ecell4/lattice/LatticeSimulator.hpp" namespace "ecell4::lattice":
    cdef cppclass Cpp_LatticeSimulator "ecell4::lattice::LatticeSimulator":
        Cpp_LatticeSimulator(
            shared_ptr[Cpp_Model], shared_ptr[Cpp_LatticeWorld]) except +
        Cpp_LatticeSimulator(
            shared_ptr[Cpp_LatticeWorld]) except +
        Integer num_steps()
        Real next_time()
        void step()
        bool step(Real& upto)
        Real t()
        Real dt()
        void set_dt(Real)
        void initialize()
        void set_alpha(Real)
        Real get_alpha()
        vector[Cpp_ReactionRule] last_reactions()
        shared_ptr[Cpp_Model] model()
        shared_ptr[Cpp_LatticeWorld] world()
        Real run(Real)
        Real run(Real, shared_ptr[Cpp_Observer])
        Real run(Real, vector[shared_ptr[Cpp_Observer]])

## LatticeSimulator
#  a python wrapper for Cpp_LatticeSimulator
cdef class LatticeSimulator:
    cdef Cpp_LatticeSimulator* thisptr

cdef LatticeSimulator LatticeSimulator_from_Cpp_LatticeSimulator(Cpp_LatticeSimulator* s)

## Cpp_LatticeFactory
#  ecell4::lattice::LatticeFactory
cdef extern from "ecell4/lattice/LatticeFactory.hpp" namespace "ecell4::lattice":
    cdef cppclass Cpp_LatticeFactory "ecell4::lattice::LatticeFactory":
        Cpp_LatticeFactory() except +
        Cpp_LatticeFactory(Real) except +
        Cpp_LatticeFactory(Real, shared_ptr[Cpp_RandomNumberGenerator]&) except +
        Cpp_LatticeWorld* create_world(string)
        Cpp_LatticeWorld* create_world(Cpp_Real3&)
        Cpp_LatticeWorld* create_world(shared_ptr[Cpp_Model])
        Cpp_LatticeSimulator* create_simulator(shared_ptr[Cpp_Model], shared_ptr[Cpp_LatticeWorld])
        Cpp_LatticeSimulator* create_simulator(shared_ptr[Cpp_LatticeWorld])

## LatticeFactory
#  a python wrapper for Cpp_LatticeFactory
cdef class LatticeFactory:
    cdef Cpp_LatticeFactory* thisptr

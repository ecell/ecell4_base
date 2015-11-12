from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *


## Cpp_ReactionInfo
cdef extern from "ecell4/spatiocyte/SpatiocyteSimulator.hpp" namespace "ecell4::spatiocyte":
    cdef cppclass Cpp_ReactionInfo "ecell4::spatiocyte::ReactionInfo":
        Cpp_ReactionInfo(Real, vector[pair[Cpp_ParticleID, Cpp_Voxel]], vector[pair[Cpp_ParticleID, Cpp_Voxel]])
        Cpp_ReactionInfo(Cpp_ReactionInfo&)
        Real t()
        vector[pair[Cpp_ParticleID, Cpp_Voxel]] reactants()
        vector[pair[Cpp_ParticleID, Cpp_Voxel]] products()

## ReactionInfo
#  a python wrapper for Cpp_ReactionInfo
cdef class ReactionInfo:
    cdef Cpp_ReactionInfo* thisptr

cdef ReactionInfo ReactionInfo_from_Cpp_ReactionInfo(Cpp_ReactionInfo* ri)

## Cpp_SpatiocyteWorld
#  ecell4::spatiocyte::SpatiocyteWorld
cdef extern from "ecell4/spatiocyte/SpatiocyteWorld.hpp" namespace "ecell4::spatiocyte":
    cdef cppclass Cpp_SpatiocyteWorld "ecell4::spatiocyte::SpatiocyteWorld":
        Cpp_SpatiocyteWorld(
            Cpp_Real3& edge_lengths, const Real& voxel_radius,
            shared_ptr[Cpp_RandomNumberGenerator] rng) except +
        Cpp_SpatiocyteWorld(
            Cpp_Real3& edge_lengths, const Real& voxel_radius) except +
        Cpp_SpatiocyteWorld(Cpp_Real3& edge_lengths) except +
        Cpp_SpatiocyteWorld(string&) except +
        Cpp_SpatiocyteWorld() except +

        void set_t(Real t)
        Real t()
        Cpp_Real3& edge_lengths()
        Real volume()
        Real voxel_volume()
        Cpp_Real3 actual_lengths()
        Real get_volume()

        pair[pair[Cpp_ParticleID, Cpp_Particle], bool] new_particle(Cpp_Particle& p)
        pair[pair[Cpp_ParticleID, Cpp_Particle], bool] new_particle(Cpp_Species& sp, Cpp_Real3& pos)
        bool remove_particle(Cpp_ParticleID& pid)
        bool remove_voxel(Cpp_ParticleID& pid)
        pair[Cpp_ParticleID, Cpp_Particle] get_particle(Cpp_ParticleID& pid)
        pair[Cpp_ParticleID, Cpp_Voxel] get_voxel(Cpp_ParticleID& pid)
        pair[Cpp_ParticleID, Cpp_Voxel] get_voxel(Integer)
        # bool on_structure(Cpp_Voxel&)
        bool on_structure(Cpp_Species&, Integer)

        void set_value(Cpp_Species&, Real)
        Real get_value(Cpp_Species&)
        Real get_value_exact(Cpp_Species&)
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
        Integer num_molecules(Cpp_Species& sp)
        Integer num_molecules_exact(Cpp_Species& sp)
        void add_molecules(Cpp_Species& sp, Integer num)
        void remove_molecules(Cpp_Species& sp, Integer num)
        # shared_ptr[Cpp_GSLRandomNumberGenerator] rng()
        Integer get_neighbor(Integer, Integer)
        Integer get_neighbor_private(Integer, Integer)
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

    cdef Cpp_SpatiocyteWorld* create_spatiocyte_world_cell_list_impl_alias(
        Cpp_Real3&, Real, Cpp_Integer3&, shared_ptr[Cpp_RandomNumberGenerator]&)
    cdef Cpp_SpatiocyteWorld* create_spatiocyte_world_vector_impl_alias(
        Cpp_Real3&, Real, shared_ptr[Cpp_RandomNumberGenerator]&)

## SpatiocyteWorld
#  a python wrapper for Cpp_SpatiocyteWorld
cdef class SpatiocyteWorld:
    cdef shared_ptr[Cpp_SpatiocyteWorld]* thisptr

cdef SpatiocyteWorld SpatiocyteWorld_from_Cpp_SpatiocyteWorld(
    shared_ptr[Cpp_SpatiocyteWorld] m)

## Cpp_SpatiocyteSimulator
#  ecell4::spatiocyte::SpatiocyteSimulator
cdef extern from "ecell4/spatiocyte/SpatiocyteSimulator.hpp" namespace "ecell4::spatiocyte":
    cdef cppclass Cpp_SpatiocyteSimulator "ecell4::spatiocyte::SpatiocyteSimulator":
        Cpp_SpatiocyteSimulator(
            shared_ptr[Cpp_Model], shared_ptr[Cpp_SpatiocyteWorld]) except +
        Cpp_SpatiocyteSimulator(
            shared_ptr[Cpp_SpatiocyteWorld]) except +
        Cpp_SpatiocyteSimulator(
            shared_ptr[Cpp_Model], shared_ptr[Cpp_SpatiocyteWorld], Real) except +
        Cpp_SpatiocyteSimulator(
            shared_ptr[Cpp_SpatiocyteWorld], Real) except +
        Integer num_steps()
        Real next_time()
        void step()
        bool step(Real& upto)
        Real t()
        void set_t(Real)
        Real dt()
        void set_dt(Real)
        void initialize()
        void set_alpha(Real)
        Real get_alpha()
        Real calculate_alpha(Cpp_ReactionRule)
        bool check_reaction()
        vector[pair[Cpp_ReactionRule, Cpp_ReactionInfo]] last_reactions()
        shared_ptr[Cpp_Model] model()
        shared_ptr[Cpp_SpatiocyteWorld] world()
        void run(Real)
        void run(Real, shared_ptr[Cpp_Observer])
        void run(Real, vector[shared_ptr[Cpp_Observer]])

## SpatiocyteSimulator
#  a python wrapper for Cpp_SpatiocyteSimulator
cdef class SpatiocyteSimulator:
    cdef Cpp_SpatiocyteSimulator* thisptr

cdef SpatiocyteSimulator SpatiocyteSimulator_from_Cpp_SpatiocyteSimulator(Cpp_SpatiocyteSimulator* s)

## Cpp_SpatiocyteFactory
#  ecell4::spatiocyte::SpatiocyteFactory
cdef extern from "ecell4/spatiocyte/SpatiocyteFactory.hpp" namespace "ecell4::spatiocyte":
    cdef cppclass Cpp_SpatiocyteFactory "ecell4::spatiocyte::SpatiocyteFactory":
        Cpp_SpatiocyteFactory() except +
        Cpp_SpatiocyteFactory(Real) except +
        Cpp_SpatiocyteFactory(Real, shared_ptr[Cpp_RandomNumberGenerator]&) except +
        Cpp_SpatiocyteFactory(Real, Real) except +
        Cpp_SpatiocyteFactory(Real, Real, shared_ptr[Cpp_RandomNumberGenerator]&) except +
        Cpp_SpatiocyteWorld* create_world()
        Cpp_SpatiocyteWorld* create_world(string)
        Cpp_SpatiocyteWorld* create_world(Cpp_Real3&)
        Cpp_SpatiocyteWorld* create_world(shared_ptr[Cpp_Model])
        Cpp_SpatiocyteSimulator* create_simulator(shared_ptr[Cpp_Model], shared_ptr[Cpp_SpatiocyteWorld])
        Cpp_SpatiocyteSimulator* create_simulator(shared_ptr[Cpp_SpatiocyteWorld])

## SpatiocyteFactory
#  a python wrapper for Cpp_SpatiocyteFactory
cdef class SpatiocyteFactory:
    cdef Cpp_SpatiocyteFactory* thisptr

from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.optional cimport optional
from ecell4.core cimport *


## CppReactionInfo
cdef extern from "ecell4/spatiocyte/SpatiocyteReactions.hpp" namespace "ecell4::spatiocyte":
    cdef cppclass CppReactionInfoItem "ecell4::spatiocyte::ReactionInfo::Item":
        CppReactionInfoItem(Cpp_ParticleID, Cpp_Species, CppVoxel)
        Cpp_ParticleID pid
        Cpp_Species species
        CppVoxel voxel

    cdef cppclass CppReactionInfo "ecell4::spatiocyte::ReactionInfo":
        CppReactionInfo(Real, vector[CppReactionInfoItem], vector[CppReactionInfoItem])
        CppReactionInfo(CppReactionInfo&)
        Real t()
        vector[CppReactionInfoItem] reactants()
        vector[CppReactionInfoItem] products()

## ReactionInfo
#  a python wrapper for CppReactionInfo
cdef class ReactionInfo:
    cdef CppReactionInfo* thisptr

cdef ReactionInfo ReactionInfo_from_Cpp_ReactionInfo(CppReactionInfo* ri)

cdef class ReactionInfoItem:
    cdef CppReactionInfoItem* thisptr

cdef ReactionInfoItem wrap_reaction_info_item(CppReactionInfoItem item)

## Cpp_SpatiocyteWorld
#  ecell4::spatiocyte::SpatiocyteWorld
cdef extern from "ecell4/spatiocyte/SpatiocyteWorld.hpp" namespace "ecell4::spatiocyte":

    cdef cppclass CppVoxel "ecell4::spatiocyte::Voxel":
        CppVoxel(CppVoxel& voxel)
        Integer num_neighbors()
        CppVoxel get_neighbor(Integer)
        Cpp_Real3 position()

    cdef cppclass Cpp_SpatiocyteWorld "ecell4::spatiocyte::SpatiocyteWorld":
        Cpp_SpatiocyteWorld(
            Cpp_Real3& edge_lengths, const Real& voxel_radius,
            shared_ptr[Cpp_RandomNumberGenerator] rng) except +
        Cpp_SpatiocyteWorld(
            Cpp_Real3& edge_lengths, const Real& voxel_radius) except +
        Cpp_SpatiocyteWorld(Cpp_Real3& edge_lengths) except +
        Cpp_SpatiocyteWorld(string&) except +
        Cpp_SpatiocyteWorld() except +

        bool has_species(Cpp_Species)

        void set_t(Real t)
        Real t()
        Cpp_Real3& edge_lengths()
        Real volume()
        Real actual_volume()
        Real voxel_volume()
        Cpp_Real3 actual_lengths()
        Real get_volume(Cpp_Species)

        optional[Cpp_ParticleID] new_particle(Cpp_Particle& p)
        optional[Cpp_ParticleID] new_particle(Cpp_Species& sp, Cpp_Real3& pos)
        bool remove_particle(Cpp_ParticleID& pid)
        bool remove_voxel(Cpp_ParticleID& pid)
        pair[Cpp_ParticleID, Cpp_Particle] get_particle(Cpp_ParticleID& pid)
        optional[Cpp_ParticleVoxel] find_voxel(Cpp_ParticleID& pid)
        pair[Cpp_ParticleID, Cpp_Species] get_voxel_at(CppVoxel)

        void set_value(Cpp_Species&, Real)
        Real get_value(Cpp_Species&)
        Real get_value_exact(Cpp_Species&)
        vector[Cpp_Species] list_species()
        vector[Cpp_Species] list_structure_species()
        vector[Cpp_Species] list_non_structure_species()
        Integer num_particles()
        Integer num_particles(Cpp_Species& sp)
        Integer num_particles_exact(Cpp_Species& sp)
        Integer num_voxels()
        Integer num_voxels(Cpp_Species& sp)
        Integer num_voxels_exact(Cpp_Species& sp)
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles()
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles(Cpp_Species& sp)
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles_exact(Cpp_Species& sp)
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_structure_particles()
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_non_structure_particles()
        bool has_particle(Cpp_ParticleID& pid)
        bool update_particle(Cpp_ParticleID& pid, Cpp_Particle& p)
        Integer num_molecules(Cpp_Species& sp)
        Integer num_molecules_exact(Cpp_Species& sp)
        void add_molecules(Cpp_Species& sp, Integer num)
        void remove_molecules(Cpp_Species& sp, Integer num)
        void save(string filename) except +
        void load(string filename)
        optional[Cpp_ParticleID] new_voxel(Cpp_Species& sp, CppVoxel pos)
        optional[Cpp_ParticleID] new_voxel_structure(Cpp_Species& sp, CppVoxel pos)
        vector[pair[Cpp_ParticleID, Cpp_ParticleVoxel]] list_voxels()
        vector[pair[Cpp_ParticleID, Cpp_ParticleVoxel]] list_voxels(Cpp_Species& sp)
        vector[pair[Cpp_ParticleID, Cpp_ParticleVoxel]] list_voxels_exact(Cpp_Species& sp)
        bool update_voxel(Cpp_ParticleID, Cpp_ParticleVoxel)
        bool has_voxel(Cpp_ParticleID)
        Real voxel_radius()

        Integer size()
        Cpp_Integer3 shape()

        void bind_to(shared_ptr[Cpp_Model])
        shared_ptr[Cpp_RandomNumberGenerator] rng()

        CppVoxel position2voxel(Cpp_Real3)

        Integer add_structure(Cpp_Species&, shared_ptr[Cpp_Shape]) except +
        void add_molecules(Cpp_Species& sp, Integer num, shared_ptr[Cpp_Shape])

        @staticmethod
        Real calculate_voxel_volume(Real)
        @staticmethod
        Cpp_Real3 calculate_hcp_lengths(Real)
        @staticmethod
        Cpp_Integer3 calculate_shape(Cpp_Real3&, Real)
        @staticmethod
        Real calculate_volume(Cpp_Real3&, Real)

    cdef Cpp_SpatiocyteWorld* create_spatiocyte_world_cell_list_impl_alias(
        Cpp_Real3&, Real, Cpp_Integer3&, shared_ptr[Cpp_RandomNumberGenerator]&)
    cdef Cpp_SpatiocyteWorld* create_spatiocyte_world_vector_impl_alias(
        Cpp_Real3&, Real, shared_ptr[Cpp_RandomNumberGenerator]&)
    cdef Cpp_SpatiocyteWorld* allocate_spatiocyte_world_square_offlattice_impl(
        Real, Real, shared_ptr[Cpp_RandomNumberGenerator]&)

## Voxel
cdef class Voxel:
    cdef CppVoxel *thisptr

cdef Voxel wrap_voxel(CppVoxel voxel)

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
        Integer num_steps()
        Real next_time()
        void step() except +
        bool step(Real& upto) except +
        Real t()
        void set_t(Real)
        Real dt()
        void set_dt(Real)
        void initialize()
        bool check_reaction()
        vector[pair[Cpp_ReactionRule, CppReactionInfo]] last_reactions()
        shared_ptr[Cpp_Model] model()
        shared_ptr[Cpp_SpatiocyteWorld] world()
        void run(Real) except +
        void run(Real, shared_ptr[Cpp_Observer]) except +
        void run(Real, vector[shared_ptr[Cpp_Observer]]) except +

## SpatiocyteSimulator
#  a python wrapper for Cpp_SpatiocyteSimulator
cdef class SpatiocyteSimulator:
    cdef Cpp_SpatiocyteSimulator* thisptr

cdef SpatiocyteSimulator SpatiocyteSimulator_from_Cpp_SpatiocyteSimulator(Cpp_SpatiocyteSimulator* s)

## Cpp_SpatiocyteFactory
#  ecell4::spatiocyte::SpatiocyteFactory
cdef extern from "ecell4/spatiocyte/SpatiocyteFactory.hpp" namespace "ecell4::spatiocyte":
    cdef cppclass Cpp_SpatiocyteFactory "ecell4::spatiocyte::SpatiocyteFactory":
        Cpp_SpatiocyteFactory(Real) except +
        Cpp_SpatiocyteWorld* create_world()
        Cpp_SpatiocyteWorld* create_world(string)
        Cpp_SpatiocyteWorld* create_world(Cpp_Real3&)
        Cpp_SpatiocyteWorld* create_world(shared_ptr[Cpp_Model])
        Cpp_SpatiocyteSimulator* create_simulator(shared_ptr[Cpp_Model], shared_ptr[Cpp_SpatiocyteWorld])
        Cpp_SpatiocyteSimulator* create_simulator(shared_ptr[Cpp_SpatiocyteWorld])
        Cpp_SpatiocyteFactory* rng_ptr(shared_ptr[Cpp_RandomNumberGenerator]&)
        @staticmethod
        Real default_voxel_radius()

## SpatiocyteFactory
#  a python wrapper for Cpp_SpatiocyteFactory
cdef class SpatiocyteFactory:
    cdef Cpp_SpatiocyteFactory* thisptr

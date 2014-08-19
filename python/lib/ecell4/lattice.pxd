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
            Cpp_Position3& edge_lengths, const Real& voxel_radius,
            shared_ptr[Cpp_RandomNumberGenerator] rng) except +
        Cpp_LatticeWorld(
            Cpp_Position3& edge_lengths, const Real& voxel_radius) except +
        Cpp_LatticeWorld(Cpp_Position3& edge_lengths) except +

        void set_t(Real t)
        Real t()
        Cpp_Position3 edge_lengths()
        Real volume()

        pair[pair[Cpp_ParticleID, Cpp_Particle], bool] new_particle(Cpp_Particle& p)
        pair[pair[Cpp_ParticleID, Cpp_Particle], bool] new_particle(Cpp_Species& sp, Cpp_Position3& pos)
        bool remove_particle(Cpp_ParticleID& pid)
        bool remove_voxel(Cpp_ParticleID& pid)
        pair[Cpp_ParticleID, Cpp_Particle] get_particle(Cpp_ParticleID& pid)
        pair[Cpp_ParticleID, Cpp_Voxel] get_voxel(Cpp_ParticleID& pid)

        Integer num_particles()
        Integer num_particles(Cpp_Species& sp)
        Integer num_voxels()
        Integer num_voxels(Cpp_Species& sp)
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles()
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles(Cpp_Species& sp)
        bool has_particle(Cpp_ParticleID& pid)
        bool update_particle(Cpp_ParticleID& pid, Cpp_Particle& p)
        # vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]] list_particles_within_radius(Cpp_Position3& pos, Real& radius)
        # vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]] list_particles_within_radius(Cpp_Position3& pos, Real& radius, Cpp_ParticleID& ignore)
        # vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]] list_particles_within_radius(Cpp_Position3& pos, Real& radius, Cpp_ParticleID& ignore1, Cpp_ParticleID& ignore2)
        # Cpp_Position3 periodic_transpose(Cpp_Position3& pos1, Cpp_Position3& pos2)
        # Cpp_Position3 apply_boundary(Cpp_Position3& pos)
        # Real distance_sq(Cpp_Position3& pos1, Cpp_Position3& pos2)
        # Real distance(Cpp_Position3& pos1, Cpp_Position3& pos2)
        # Real volume()
        # # bool has_species(Cpp_Species& sp)
        Integer num_molecules()
        Integer num_molecules(Cpp_Species& sp)
        void add_molecules(Cpp_Species& sp, Integer num)
        void remove_molecules(Cpp_Species& sp, Integer num)
        # shared_ptr[Cpp_GSLRandomNumberGenerator] rng()
        void save(string filename)
        void load(string filename)
        pair[pair[Cpp_ParticleID, Cpp_Voxel], bool] new_voxel(Cpp_Voxel& p)
        pair[pair[Cpp_ParticleID, Cpp_Voxel], bool] new_voxel(Cpp_Species& sp, Integer pos)
        vector[pair[Cpp_ParticleID, Cpp_Voxel]] list_voxels(Cpp_Species& sp)
        bool update_voxel(Cpp_ParticleID, Cpp_Voxel)
        bool has_voxel(Cpp_ParticleID)
        Real voxel_radius()
        Integer col_size()
        Integer row_size()
        Integer layer_size()
        Integer size()
        void bind_to(shared_ptr[Cpp_Model])
        Cpp_Position3 coordinate2position(Integer)
        Integer position2coordinate(Cpp_Position3)
        shared_ptr[Cpp_RandomNumberGenerator] rng()

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
        Integer num_steps()
        Real next_time()
        void step()
        bool step(Real& upto)
        Real t()
        Real dt()
        void set_dt(Real)
        void initialize()
        vector[Cpp_ReactionRule] last_reactions()
        shared_ptr[Cpp_Model] model()
        shared_ptr[Cpp_LatticeWorld] world()
        void run(Real)
        void run(Real, vector[shared_ptr[Cpp_Observer]])

## LatticeSimulator
#  a python wrapper for Cpp_LatticeSimulator
cdef class LatticeSimulator:
    cdef Cpp_LatticeSimulator* thisptr


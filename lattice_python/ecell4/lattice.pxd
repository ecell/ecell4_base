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
            shared_ptr[Cpp_GSLRandomNumberGenerator] rng) except +
        Cpp_LatticeWorld(
            Cpp_Position3& edge_lengths, const Real& voxel_radius) except +
        Cpp_LatticeWorld(Cpp_Position3& edge_lengths) except +

        void set_t(Real t)
        Real t()
        Cpp_Position3 edge_lengths()
        Real volume()

        # Cpp_ParticleID new_particle(Cpp_Particle& p)
        # Integer num_particles()
        Integer num_particles(Cpp_Species& sp)
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles()
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles(Cpp_Species& sp)
        bool has_particle(Cpp_ParticleID& pid)
        bool update_particle(Cpp_ParticleID& pid, Cpp_Particle& p)
        # pair[Cpp_ParticleID, Cpp_Particle] get_particle(Cpp_ParticleID& pid)
        # void remove_particle(Cpp_ParticleID& pid)
        # vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]] list_particles_within_radius(Cpp_Position3& pos, Real& radius)
        # vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]] list_particles_within_radius(Cpp_Position3& pos, Real& radius, Cpp_ParticleID& ignore)
        # vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]] list_particles_within_radius(Cpp_Position3& pos, Real& radius, Cpp_ParticleID& ignore1, Cpp_ParticleID& ignore2)
        # Cpp_Position3 periodic_transpose(Cpp_Position3& pos1, Cpp_Position3& pos2)
        # Cpp_Position3 apply_boundary(Cpp_Position3& pos)
        # Real distance_sq(Cpp_Position3& pos1, Cpp_Position3& pos2)
        # Real distance(Cpp_Position3& pos1, Cpp_Position3& pos2)
        # Real volume()
        # # bool has_species(Cpp_Species& sp)
        Integer num_molecules(Cpp_Species& sp)
        void add_molecules(Cpp_Species& sp, Integer num)
        void remove_molecules(Cpp_Species& sp, Integer num)
        # shared_ptr[Cpp_GSLRandomNumberGenerator] rng()
        void save(string filename)
        void load(string filename)
        vector[pair[Cpp_ParticleID, Cpp_Voxel]] list_voxels(Cpp_Species& sp)
        Real voxel_radius()
        Integer col_size()
        Integer row_size()
        Integer layer_size()
        Integer size()
        void bind_to(shared_ptr[Cpp_NetworkModel])

## LatticeWorld
#  a python wrapper for Cpp_LatticeWorld
cdef class LatticeWorld:
    cdef shared_ptr[Cpp_LatticeWorld]* thisptr

## Cpp_LatticeSimulator
#  ecell4::lattice::LatticeSimulator
cdef extern from "ecell4/lattice/LatticeSimulator.hpp" namespace "ecell4::lattice":
    cdef cppclass Cpp_LatticeSimulator "ecell4::lattice::LatticeSimulator":
        Cpp_LatticeSimulator(
            shared_ptr[Cpp_NetworkModel], shared_ptr[Cpp_LatticeWorld]) except +
        Integer num_steps()
        void step()
        bool step(Real& upto)
        Real t()
        Real dt()
        # void initialize()

## LatticeSimulator
#  a python wrapper for Cpp_LatticeSimulator
cdef class LatticeSimulator:
    cdef Cpp_LatticeSimulator* thisptr


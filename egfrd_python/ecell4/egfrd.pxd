from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *


## Cpp_EGFRDWorld
#  ecell4::egfrd::EGFRDWorld
cdef extern from "ecell4/egfrd/EGFRDWorld.hpp" namespace "ecell4::egfrd":
    cdef cppclass Cpp_EGFRDWorld "ecell4::egfrd::EGFRDWorld":
        Cpp_EGFRDWorld(
            Real world_size, Integer matrix_size,
            shared_ptr[Cpp_GSLRandomNumberGenerator] rng) except +
        Cpp_ParticleID new_particle(Cpp_Particle& p)
        void set_t(Real t)
        Real t()
        Cpp_Position3 edge_lengths()
        Integer num_particles()
        Integer num_particles(Cpp_Species& sp)
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles()
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles(Cpp_Species& sp)
        bool has_particle(Cpp_ParticleID& pid)
        bool update_particle(Cpp_ParticleID& pid, Cpp_Particle& p)
        pair[Cpp_ParticleID, Cpp_Particle] get_particle(Cpp_ParticleID& pid)
        void remove_particle(Cpp_ParticleID& pid)
        vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]] list_particles_within_radius(Cpp_Position3& pos, Real& radius)
        vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]] list_particles_within_radius(Cpp_Position3& pos, Real& radius, Cpp_ParticleID& ignore)
        vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]] list_particles_within_radius(Cpp_Position3& pos, Real& radius, Cpp_ParticleID& ignore1, Cpp_ParticleID& ignore2)
        Cpp_Position3 periodic_transpose(Cpp_Position3& pos1, Cpp_Position3& pos2)
        Cpp_Position3 apply_boundary(Cpp_Position3& pos)
        Real distance_sq(Cpp_Position3& pos1, Cpp_Position3& pos2)
        Real distance(Cpp_Position3& pos1, Cpp_Position3& pos2)
        Real volume()
        bool has_species(Cpp_Species& sp)
        Integer num_molecules(Cpp_Species& sp)
        void add_molecules(Cpp_Species& sp, Integer num)
        shared_ptr[GSLRandomNumberGenerator] rng()
        void save(string filename)

## EGFRDWorld
#  a python wrapper for Cpp_EGFRDWorld
cdef class EGFRDWorld:
    cdef shared_ptr[Cpp_EGFRDWorld]* thisptr

## Cpp_EGFRDSimulator
#  ecell4::egfrd::EGFRDSimulator
cdef extern from "ecell4/egfrd/EGFRDSimulatorWrapper.hpp" namespace "ecell4::egfrd":
    cdef cppclass Cpp_EGFRDSimulatorWrapper "ecell4::egfrd::EGFRDSimulatorWrapper":
        # Cpp_EGFRDSimulatorWrapper(
        #     shared_ptr[Cpp_NetworkModel], shared_ptr[Cpp_EGFRDWorld],
        #     Integer dissociation_retry_moves) except +
        Cpp_EGFRDSimulatorWrapper(
            shared_ptr[Cpp_NetworkModel], shared_ptr[Cpp_EGFRDWorld]) except +
        Integer num_steps()
        void step()
        bool step(Real& upto)
        Real t()
        Real dt()
        void initialize()

## EGFRDSimulator
#  a python wrapper for Cpp_EGFRDSimulatorWrapper
cdef class EGFRDSimulatorWrapper:
    cdef Cpp_EGFRDSimulatorWrapper* thisptr

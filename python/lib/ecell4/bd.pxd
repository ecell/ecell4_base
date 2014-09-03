from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *


## Cpp_BDWorld
#  ecell4::bd::BDWorld
cdef extern from "ecell4/bd/BDWorld.hpp" namespace "ecell4::bd":
    cdef cppclass Cpp_BDWorld "ecell4::bd::BDWorld":
        Cpp_BDWorld() except +
        Cpp_BDWorld(string& edge_lengths) except +
        Cpp_BDWorld(Cpp_Position3& edge_lengths) except +
        Cpp_BDWorld(
            Cpp_Position3& edge_lengths,
            shared_ptr[Cpp_RandomNumberGenerator] rng) except +
        pair[pair[Cpp_ParticleID, Cpp_Particle], bool] new_particle(Cpp_Particle& p)
        pair[pair[Cpp_ParticleID, Cpp_Particle], bool] new_particle(Cpp_Species& sp, Cpp_Position3& pos)
        void set_t(Real t)
        Real t()
        Cpp_Position3 edge_lengths()
        Integer num_particles()
        Integer num_particles(Cpp_Species& sp)
        Integer num_particles_exact(Cpp_Species& sp)
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles()
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles(Cpp_Species& sp)
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles_exact(Cpp_Species& sp)
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
        # bool has_species(Cpp_Species& sp)
        Integer num_molecules(Cpp_Species& sp)
        Integer num_molecules_exact(Cpp_Species& sp)
        void add_molecules(Cpp_Species& sp, Integer num)
        void remove_molecules(Cpp_Species& sp, Integer num)
        void save(string filename)
        void load(string filename)
        void bind_to(shared_ptr[Cpp_Model])
        shared_ptr[Cpp_RandomNumberGenerator] rng()

## BDWorld
#  a python wrapper for Cpp_BDWorld
cdef class BDWorld:
    cdef shared_ptr[Cpp_BDWorld]* thisptr

cdef BDWorld BDWorld_from_Cpp_BDWorld(shared_ptr[Cpp_BDWorld] m)

## Cpp_BDSimulator
#  ecell4::bd::BDSimulator
cdef extern from "ecell4/bd/BDSimulator.hpp" namespace "ecell4::bd":
    cdef cppclass Cpp_BDSimulator "ecell4::bd::BDSimulator":
        # Cpp_BDSimulator(
        #     shared_ptr[Cpp_NetworkModel], shared_ptr[Cpp_BDWorld],
        #     Integer dissociation_retry_moves) except +
        Cpp_BDSimulator(
            shared_ptr[Cpp_Model], shared_ptr[Cpp_BDWorld]) except +
        Integer num_steps()
        void step()
        bool step(Real& upto)
        Real t()
        Real dt()
        void set_dt(Real& dt)
        vector[Cpp_ReactionRule] last_reactions()
        Real next_time()
        void initialize()
        shared_ptr[Cpp_Model] model()
        shared_ptr[Cpp_BDWorld] world()
        void run(Real)
        void run(Real, vector[shared_ptr[Cpp_Observer]])

## BDSimulator
#  a python wrapper for Cpp_BDSimulator
cdef class BDSimulator:
    cdef Cpp_BDSimulator* thisptr

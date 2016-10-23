from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *


## Cpp_ReactionInfo
cdef extern from "ecell4/bd/BDPropagator.hpp" namespace "ecell4::bd":
    cdef cppclass Cpp_ReactionInfo "ecell4::bd::ReactionInfo":
        Cpp_ReactionInfo(Real, vector[pair[Cpp_ParticleID, Cpp_Particle]], vector[pair[Cpp_ParticleID, Cpp_Particle]])
        Cpp_ReactionInfo(Cpp_ReactionInfo&)
        Real t()
        vector[pair[Cpp_ParticleID, Cpp_Particle]] reactants()
        vector[pair[Cpp_ParticleID, Cpp_Particle]] products()

## ReactionInfo
#  a python wrapper for Cpp_ReactionInfo
cdef class ReactionInfo:
    cdef Cpp_ReactionInfo* thisptr

cdef ReactionInfo ReactionInfo_from_Cpp_ReactionInfo(Cpp_ReactionInfo* ri)

## Cpp_BDWorld
#  ecell4::bd::BDWorld
cdef extern from "ecell4/bd/BDWorld.hpp" namespace "ecell4::bd":
    cdef cppclass Cpp_BDWorld "ecell4::bd::BDWorld":
        Cpp_BDWorld() except +
        Cpp_BDWorld(string& filename) except +
        Cpp_BDWorld(Cpp_Real3& edge_lengths) except +
        Cpp_BDWorld(
            Cpp_Real3& edge_lengths,
            Cpp_Integer3& matrix_sizes) except +
        Cpp_BDWorld(
            Cpp_Real3& edge_lengths,
            Cpp_Integer3& matrix_sizes,
            shared_ptr[Cpp_RandomNumberGenerator] rng) except +

        pair[pair[Cpp_ParticleID, Cpp_Particle], bool] new_particle(Cpp_Particle& p)
        pair[pair[Cpp_ParticleID, Cpp_Particle], bool] new_particle(Cpp_Species& sp, Cpp_Real3& pos)
        void set_t(Real t)
        Real t()
        Cpp_Real3& edge_lengths()
        Cpp_Real3 actual_lengths()
        Real get_value(Cpp_Species&)
        Real get_value_exact(Cpp_Species&)
        vector[Cpp_Species] list_species()
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
        vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]] list_particles_within_radius(Cpp_Real3& pos, Real& radius)
        vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]] list_particles_within_radius(Cpp_Real3& pos, Real& radius, Cpp_ParticleID& ignore)
        vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]] list_particles_within_radius(Cpp_Real3& pos, Real& radius, Cpp_ParticleID& ignore1, Cpp_ParticleID& ignore2)
        Cpp_Real3 periodic_transpose(Cpp_Real3& pos1, Cpp_Real3& pos2)
        Cpp_Real3 apply_boundary(Cpp_Real3& pos)
        Real distance_sq(Cpp_Real3& pos1, Cpp_Real3& pos2)
        Real distance(Cpp_Real3& pos1, Cpp_Real3& pos2)
        Real volume()
        # bool has_species(Cpp_Species& sp)
        Integer num_molecules(Cpp_Species& sp)
        Integer num_molecules_exact(Cpp_Species& sp)
        void add_molecules(Cpp_Species& sp, Integer num)
        void add_molecules(Cpp_Species& sp, Integer num, shared_ptr[Cpp_Shape])
        void remove_molecules(Cpp_Species& sp, Integer num)
        void save(string filename) except +
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
        Cpp_BDSimulator(
            shared_ptr[Cpp_BDWorld]) except +
        Cpp_BDSimulator(
            shared_ptr[Cpp_Model], shared_ptr[Cpp_BDWorld], Real) except +
        Cpp_BDSimulator(
            shared_ptr[Cpp_BDWorld], Real) except +
        Integer num_steps()
        void step() except +
        bool step(Real& upto) except +
        Real t()
        Real dt()
        void set_t(Real& t)
        void set_dt(Real& dt)
        vector[pair[Cpp_ReactionRule, Cpp_ReactionInfo]] last_reactions()
        bool check_reaction()
        Real next_time()
        void initialize()
        shared_ptr[Cpp_Model] model()
        shared_ptr[Cpp_BDWorld] world()
        void run(Real) except +
        void run(Real, shared_ptr[Cpp_Observer]) except +
        void run(Real, vector[shared_ptr[Cpp_Observer]]) except +

## BDSimulator
#  a python wrapper for Cpp_BDSimulator
cdef class BDSimulator:
    cdef Cpp_BDSimulator* thisptr

cdef BDSimulator BDSimulator_from_Cpp_BDSimulator(Cpp_BDSimulator* s)

## Cpp_BDFactory
#  ecell4::bd::BDFactory
cdef extern from "ecell4/bd/BDFactory.hpp" namespace "ecell4::bd":
    cdef cppclass Cpp_BDFactory "ecell4::bd::BDFactory":
        Cpp_BDFactory(Cpp_Integer3&, Real) except +
        Cpp_BDWorld* create_world()
        Cpp_BDWorld* create_world(string)
        Cpp_BDWorld* create_world(Cpp_Real3&)
        Cpp_BDWorld* create_world(shared_ptr[Cpp_Model])
        Cpp_BDSimulator* create_simulator(shared_ptr[Cpp_Model], shared_ptr[Cpp_BDWorld])
        Cpp_BDSimulator* create_simulator(shared_ptr[Cpp_BDWorld])
        Cpp_BDFactory* rng_ptr(shared_ptr[Cpp_RandomNumberGenerator]&)
        @staticmethod
        Cpp_Integer3 default_matrix_sizes()
        @staticmethod
        Real default_bd_dt_factor()

## BDFactory
#  a python wrapper for Cpp_BDFactory
cdef class BDFactory:
    cdef Cpp_BDFactory* thisptr

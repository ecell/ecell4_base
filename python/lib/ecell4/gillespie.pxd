from libcpp.string cimport string
from libcpp cimport bool

from core cimport *


## Cpp_ReactionInfo
cdef extern from "ecell4/gillespie/GillespieSimulator.hpp" namespace "ecell4::gillespie":
    cdef cppclass Cpp_ReactionInfo "ecell4::gillespie::ReactionInfo":
        Cpp_ReactionInfo(Real, vector[Cpp_Species], vector[Cpp_Species])
        Cpp_ReactionInfo(Cpp_ReactionInfo&)
        Real t()
        vector[Cpp_Species] reactants()
        vector[Cpp_Species] products()

## ReactionInfo
#  a python wrapper for Cpp_ReactionInfo
cdef class ReactionInfo:
    cdef Cpp_ReactionInfo* thisptr

cdef ReactionInfo ReactionInfo_from_Cpp_ReactionInfo(Cpp_ReactionInfo* ri)

## Cpp_GillespieWorld
#  ecell4::gillespie::GillespieWorld
cdef extern from "ecell4/gillespie/GillespieWorld.hpp" namespace "ecell4::gillespie":
    cdef cppclass Cpp_GillespieWorld "ecell4::gillespie::GillespieWorld":
        Cpp_GillespieWorld() except +
        Cpp_GillespieWorld(Cpp_Real3&) except +
        Cpp_GillespieWorld(string&) except +
        Cpp_GillespieWorld(Cpp_Real3&, shared_ptr[Cpp_RandomNumberGenerator]) except +
        void set_t(Real)
        Real t()
        Real volume()
        void reset(Cpp_Real3&)
        Cpp_Real3 edge_lengths()
        Integer num_molecules(Cpp_Species &)
        Integer num_molecules_exact(Cpp_Species &)
        vector[Cpp_Species] list_species()
        void add_molecules(Cpp_Species &sp, Integer &num)
        void add_molecules(Cpp_Species &sp, Integer &num, shared_ptr[Cpp_Shape])
        void remove_molecules(Cpp_Species &sp, Integer &num)
        void save(string) except +
        void load(string)
        void bind_to(shared_ptr[Cpp_Model])
        shared_ptr[Cpp_RandomNumberGenerator] rng()
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles()
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles(Cpp_Species& sp)
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles_exact(Cpp_Species& sp)

## GillespieWorld
#  a python wrapper for Cpp_GillespieWorld
cdef class GillespieWorld:
    cdef shared_ptr[Cpp_GillespieWorld]* thisptr

cdef GillespieWorld GillespieWorld_from_Cpp_GillespieWorld(
    shared_ptr[Cpp_GillespieWorld] m)

## Cpp_GillespieSimulator
#  ecell4::gillespie::GillespieSimulator
cdef extern from "ecell4/gillespie/GillespieSimulator.hpp" namespace "ecell4::gillespie":
    cdef cppclass Cpp_GillespieSimulator "ecell4::gillespie::GillespieSimulator":
        Cpp_GillespieSimulator(
            shared_ptr[Cpp_Model], shared_ptr[Cpp_GillespieWorld]) except +
        Cpp_GillespieSimulator(
            shared_ptr[Cpp_GillespieWorld]) except +
        Integer num_steps()
        void step()
        bool step(Real)
        Real t()
        void set_t(Real)
        void set_dt(Real)
        Real dt()
        Real next_time()
        bool check_reaction()
        vector[pair[Cpp_ReactionRule, Cpp_ReactionInfo]] last_reactions()
        void initialize()
        # Cpp_GSLRandomNumberGenerator& rng()
        shared_ptr[Cpp_Model] model()
        shared_ptr[Cpp_GillespieWorld] world()
        void run(Real)
        void run(Real, shared_ptr[Cpp_Observer])
        void run(Real, vector[shared_ptr[Cpp_Observer]])

## GillespieSimulator
#  a python wrapper for Cpp_GillespieSimulator
cdef class GillespieSimulator:
    cdef Cpp_GillespieSimulator* thisptr

cdef GillespieSimulator GillespieSimulator_from_Cpp_GillespieSimulator(Cpp_GillespieSimulator* s)

## Cpp_GillespieFactory
#  ecell4::gillespie::GillespieFactory
cdef extern from "ecell4/gillespie/GillespieFactory.hpp" namespace "ecell4::gillespie":
    cdef cppclass Cpp_GillespieFactory "ecell4::gillespie::GillespieFactory":
        Cpp_GillespieFactory() except +
        Cpp_GillespieFactory(shared_ptr[Cpp_RandomNumberGenerator]) except +
        Cpp_GillespieWorld* create_world(string)
        Cpp_GillespieWorld* create_world(Cpp_Real3&)
        Cpp_GillespieWorld* create_world(shared_ptr[Cpp_Model])
        Cpp_GillespieSimulator* create_simulator(shared_ptr[Cpp_Model], shared_ptr[Cpp_GillespieWorld])
        Cpp_GillespieSimulator* create_simulator(shared_ptr[Cpp_GillespieWorld])

## GillespieFactory
#  a python wrapper for Cpp_GillespieFactory
cdef class GillespieFactory:
    cdef Cpp_GillespieFactory* thisptr

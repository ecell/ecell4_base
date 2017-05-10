from libcpp.string cimport string
from libcpp cimport bool

from core cimport *


## Cpp_ReactionInfo
cdef extern from "ecell4/egfrd/egfrd.hpp" namespace "ecell4::egfrd":
    cdef cppclass Cpp_ReactionInfo "ecell4::egfrd::ReactionInfo":
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

## Cpp_EGFRDWorld
#  ecell4::egfrd::EGFRDWorld
cdef extern from "ecell4/egfrd/egfrd.hpp" namespace "ecell4::egfrd":
    cdef cppclass Cpp_EGFRDWorld "ecell4::egfrd::EGFRDWorld":
        Cpp_EGFRDWorld() except +
        Cpp_EGFRDWorld(Cpp_Real3&) except +
        Cpp_EGFRDWorld(Cpp_Real3&, Cpp_Integer3&) except +
        Cpp_EGFRDWorld(
            Cpp_Real3&, Cpp_Integer3&,
            shared_ptr[Cpp_RandomNumberGenerator]&) except +
        #     shared_ptr[Cpp_GSLRandomNumberGenerator]&) except +
        Cpp_EGFRDWorld(string&) except +
        pair[pair[Cpp_ParticleID, Cpp_Particle], bool] new_particle(Cpp_Particle& p)
        pair[pair[Cpp_ParticleID, Cpp_Particle], bool] new_particle(Cpp_Species& sp, Cpp_Real3& pos)
        void set_t(Real t)
        Real t()
        Cpp_Real3& edge_lengths()
        Cpp_Real3 actual_lengths()
        void set_value(Cpp_Species&, Real)
        Real get_value(Cpp_Species&)
        Real get_value_exact(Cpp_Species&)
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
        # Cpp_Real3 periodic_transpose(Cpp_Real3& pos1, Cpp_Real3& pos2)
        Cpp_Real3 apply_boundary(Cpp_Real3& pos)
        # Real distance_sq(Cpp_Real3& pos1, Cpp_Real3& pos2)
        Real distance(Cpp_Real3& pos1, Cpp_Real3& pos2)
        Real volume()
        bool has_species(Cpp_Species& sp)
        vector[Cpp_Species] list_species()
        Integer num_molecules(Cpp_Species& sp)
        Integer num_molecules_exact(Cpp_Species& sp)
        void add_molecules(Cpp_Species& sp, Integer num)
        void add_molecules(Cpp_Species& sp, Integer num, shared_ptr[Cpp_Shape])
        void remove_molecules(Cpp_Species& sp, Integer num)
        void save(string filename) except +
        void load(string filename) except +
        void bind_to(shared_ptr[Cpp_Model])
        shared_ptr[Cpp_RandomNumberGenerator] rng()

    cdef cppclass Cpp_EGFRDSimulator "ecell4::egfrd::EGFRDSimulator":
        #XXX: be carefull about the order of arguments
        Cpp_EGFRDSimulator(
            shared_ptr[Cpp_EGFRDWorld]&, shared_ptr[Cpp_Model]&) except +
        Cpp_EGFRDSimulator(
            shared_ptr[Cpp_EGFRDWorld]&, shared_ptr[Cpp_Model]&,
            Real) except +
        Cpp_EGFRDSimulator(
            shared_ptr[Cpp_EGFRDWorld]&, shared_ptr[Cpp_Model]&,
            Real, Integer) except +
        Cpp_EGFRDSimulator(
            shared_ptr[Cpp_EGFRDWorld]&, shared_ptr[Cpp_Model]&,
            Real, Integer, Real) except +
        Cpp_EGFRDSimulator(
            shared_ptr[Cpp_EGFRDWorld]&) except +
        Cpp_EGFRDSimulator(
            shared_ptr[Cpp_EGFRDWorld]&, Real) except +
        Cpp_EGFRDSimulator(
            shared_ptr[Cpp_EGFRDWorld]&, Real, Integer) except +
        Cpp_EGFRDSimulator(
            shared_ptr[Cpp_EGFRDWorld]&, Real, Integer, Real) except +
        Integer num_steps()
        void step() except +
        bool step(Real) except +
        Real t()
        void set_t(Real)
        void set_dt(Real)
        Real dt()
        Real next_time()
        vector[pair[Cpp_ReactionRule, Cpp_ReactionInfo]] last_reactions()
        bool check_reaction()
        void initialize()
        # Cpp_GSLRandomNumberGenerator& rng()
        shared_ptr[Cpp_Model] model()
        shared_ptr[Cpp_EGFRDWorld] world()
        void run(Real) except +
        void run(Real, shared_ptr[Cpp_Observer]) except +
        void run(Real, vector[shared_ptr[Cpp_Observer]]) except +
        void set_paranoiac(bool)

    cdef cppclass Cpp_EGFRDFactory "ecell4::egfrd::EGFRDFactory":
        Cpp_EGFRDFactory(Cpp_Integer3&, Real, Integer, Real) except +
        Cpp_EGFRDWorld* create_world()
        Cpp_EGFRDWorld* create_world(string)
        Cpp_EGFRDWorld* create_world(Cpp_Real3&)
        Cpp_EGFRDWorld* create_world(shared_ptr[Cpp_Model])
        Cpp_EGFRDSimulator* create_simulator(shared_ptr[Cpp_Model], shared_ptr[Cpp_EGFRDWorld])
        Cpp_EGFRDSimulator* create_simulator(shared_ptr[Cpp_EGFRDWorld])
        Cpp_EGFRDFactory* rng_ptr(shared_ptr[Cpp_RandomNumberGenerator]&)
        @staticmethod
        Cpp_Integer3 default_matrix_sizes()
        @staticmethod
        Real default_bd_dt_factor()
        @staticmethod
        Integer default_dissociation_retry_moves()
        @staticmethod
        Real default_user_max_shell_size()

    cdef cppclass Cpp_BDFactory "ecell4::egfrd::BDFactory":
        Cpp_BDFactory(Cpp_Integer3&, Real, Integer) except +
        Cpp_EGFRDWorld* create_world()
        Cpp_EGFRDWorld* create_world(string)
        Cpp_EGFRDWorld* create_world(Cpp_Real3&)
        Cpp_EGFRDWorld* create_world(shared_ptr[Cpp_Model])
        Cpp_BDSimulator* create_simulator(shared_ptr[Cpp_Model], shared_ptr[Cpp_EGFRDWorld])
        Cpp_BDSimulator* create_simulator(shared_ptr[Cpp_EGFRDWorld])
        Cpp_BDFactory* rng_ptr(shared_ptr[Cpp_RandomNumberGenerator]&)
        @staticmethod
        Cpp_Integer3 default_matrix_sizes()
        @staticmethod
        Real default_bd_dt_factor()
        @staticmethod
        Integer default_dissociation_retry_moves()

    cdef cppclass Cpp_BDSimulator "ecell4::egfrd::BDSimulator":
        #XXX: be carefull about the order of arguments
        Cpp_BDSimulator(
            shared_ptr[Cpp_EGFRDWorld]&, shared_ptr[Cpp_Model]&) except +
        Cpp_BDSimulator(
            shared_ptr[Cpp_EGFRDWorld]&, shared_ptr[Cpp_Model]&,
            Real) except +
        Cpp_BDSimulator(
            shared_ptr[Cpp_EGFRDWorld]&, shared_ptr[Cpp_Model]&,
            Real, Integer) except +
        Cpp_BDSimulator(shared_ptr[Cpp_EGFRDWorld]&) except +
        Cpp_BDSimulator(shared_ptr[Cpp_EGFRDWorld]&, Real) except +
        Cpp_BDSimulator(shared_ptr[Cpp_EGFRDWorld]&, Real, Integer) except +
        Integer num_steps()
        void step() except +
        bool step(Real) except +
        Real t()
        void set_t(Real)
        void set_dt(Real)
        Real dt()
        Real next_time()
        vector[pair[Cpp_ReactionRule, Cpp_ReactionInfo]] last_reactions()
        bool check_reaction()
        void initialize()
        # Cpp_GSLRandomNumberGenerator& rng()
        shared_ptr[Cpp_Model] model()
        shared_ptr[Cpp_EGFRDWorld] world()
        void run(Real) except +
        void run(Real, shared_ptr[Cpp_Observer]) except +
        void run(Real, vector[shared_ptr[Cpp_Observer]]) except +
        Real dt_factor()

        void add_potential(Cpp_Species&, shared_ptr[Cpp_Shape], Real) except +
        void add_potential(Cpp_Species&, shared_ptr[Cpp_Shape]) except +
        void add_potential(Cpp_Species&, Real) except +

cdef class EGFRDWorld:
    cdef shared_ptr[Cpp_EGFRDWorld]* thisptr

cdef class EGFRDSimulator:
    cdef Cpp_EGFRDSimulator* thisptr

cdef class EGFRDFactory:
    cdef Cpp_EGFRDFactory* thisptr

cdef class BDSimulator:
    cdef Cpp_BDSimulator* thisptr

cdef class BDFactory:
    cdef Cpp_BDFactory* thisptr

cdef EGFRDWorld EGFRDWorld_from_Cpp_EGFRDWorld(
    shared_ptr[Cpp_EGFRDWorld] m)

cdef EGFRDSimulator EGFRDSimulator_from_Cpp_EGFRDSimulator(
    Cpp_EGFRDSimulator* s)

cdef BDSimulator BDSimulator_from_Cpp_BDSimulator(
    Cpp_BDSimulator* s)

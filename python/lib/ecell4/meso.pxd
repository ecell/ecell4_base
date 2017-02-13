from libcpp.string cimport string
from libcpp cimport bool

from core cimport *


## Cpp_ReactionInfo
cdef extern from "ecell4/meso/MesoscopicSimulator.hpp" namespace "ecell4::meso":
    cdef cppclass Cpp_ReactionInfo "ecell4::meso::ReactionInfo":
        Cpp_ReactionInfo(Real, vector[Cpp_Species], vector[Cpp_Species], Integer)
        Cpp_ReactionInfo(Cpp_ReactionInfo&)
        Real t()
        vector[Cpp_Species] reactants()
        vector[Cpp_Species] products()
        Integer coordinate()

## ReactionInfo
#  a python wrapper for Cpp_ReactionInfo
cdef class ReactionInfo:
    cdef Cpp_ReactionInfo* thisptr

cdef ReactionInfo ReactionInfo_from_Cpp_ReactionInfo(Cpp_ReactionInfo* ri)

## Cpp_MesoscopicWorld
#  ecell4::meso::MesoscopicWorld
cdef extern from "ecell4/meso/MesoscopicWorld.hpp" namespace "ecell4::meso":
    cdef cppclass Cpp_MesoscopicWorld "ecell4::meso::MesoscopicWorld":
        Cpp_MesoscopicWorld() except +
        Cpp_MesoscopicWorld(Cpp_Real3&) except +
        Cpp_MesoscopicWorld(string&) except +
        Cpp_MesoscopicWorld(Cpp_Real3&, Cpp_Integer3&) except +
        Cpp_MesoscopicWorld(Cpp_Real3&, Cpp_Integer3&,
            shared_ptr[Cpp_RandomNumberGenerator]) except +
        Cpp_MesoscopicWorld(Cpp_Real3&, Real&) except +
        Cpp_MesoscopicWorld(Cpp_Real3&, Real&,
            shared_ptr[Cpp_RandomNumberGenerator]) except +
        void set_t(Real)
        Real t()
        Real volume()
        Real subvolume()
        Integer num_subvolumes()
        Integer num_subvolumes(Cpp_Species&)
        void reset(Cpp_Real3&)
        Integer global2coord(Cpp_Integer3&)
        Cpp_Integer3 coord2global(Integer)
        Cpp_Integer3 position2global(Cpp_Real3&)
        Integer position2coordinate(Cpp_Real3&)
        Cpp_Real3& edge_lengths()
        Cpp_Real3 actual_lengths()
        Cpp_Integer3 matrix_sizes()
        Cpp_Real3 subvolume_edge_lengths()
        void set_value(Cpp_Species&, Real)
        Real get_value(Cpp_Species&)
        Real get_value_exact(Cpp_Species&)
        Integer num_molecules(Cpp_Species &)
        Integer num_molecules_exact(Cpp_Species &)
        Integer num_molecules(Cpp_Species &, Integer)
        Integer num_molecules_exact(Cpp_Species &, Integer)
        Integer num_molecules(Cpp_Species &, Cpp_Integer3)
        Integer num_molecules_exact(Cpp_Species &, Cpp_Integer3)
        vector[Cpp_Species] list_species()
        vector[Integer] list_coordinates(Cpp_Species&)
        vector[Integer] list_coordinates_exact(Cpp_Species&)
        void add_molecules(Cpp_Species &sp, Integer &num, Integer)
        void remove_molecules(Cpp_Species &sp, Integer &num, Integer)
        void add_molecules(Cpp_Species &sp, Integer &num, Cpp_Integer3)
        void remove_molecules(Cpp_Species &sp, Integer &num, Cpp_Integer3)
        void add_molecules(Cpp_Species &sp, Integer &num)
        void add_molecules(Cpp_Species &sp, Integer &num, shared_ptr[Cpp_Shape])
        void add_structure(Cpp_Species&, shared_ptr[Cpp_Shape]) except+
        Real get_volume(Cpp_Species&)
        Real get_occupancy(Cpp_Species&, Integer&)
        Real get_occupancy(Cpp_Species&, Cpp_Integer3&)
        bool on_structure(Cpp_Species&, Cpp_Integer3&)
        bool check_structure(Cpp_Species&, Cpp_Integer3&)
        bool has_structure(Cpp_Species&)
        void remove_molecules(Cpp_Species &sp, Integer &num)
        void save(string) except +
        void load(string) except +
        void bind_to(shared_ptr[Cpp_Model])
        shared_ptr[Cpp_RandomNumberGenerator] rng()
        pair[pair[Cpp_ParticleID, Cpp_Particle], bool] new_particle(Cpp_Particle& p)
        pair[pair[Cpp_ParticleID, Cpp_Particle], bool] new_particle(Cpp_Species& sp, Cpp_Real3& pos)
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles()
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles(Cpp_Species& sp)
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles_exact(Cpp_Species& sp)

## MesoscopicWorld
#  a python wrapper for Cpp_MesoscopicWorld
cdef class MesoscopicWorld:
    cdef shared_ptr[Cpp_MesoscopicWorld]* thisptr

cdef MesoscopicWorld MesoscopicWorld_from_Cpp_MesoscopicWorld(
    shared_ptr[Cpp_MesoscopicWorld] m)

## Cpp_MesoscopicSimulator
#  ecell4::meso::MesoscopicSimulator
cdef extern from "ecell4/meso/MesoscopicSimulator.hpp" namespace "ecell4::meso":
    cdef cppclass Cpp_MesoscopicSimulator "ecell4::meso::MesoscopicSimulator":
        Cpp_MesoscopicSimulator(
            shared_ptr[Cpp_Model], shared_ptr[Cpp_MesoscopicWorld]) except +
        Cpp_MesoscopicSimulator(
            shared_ptr[Cpp_MesoscopicWorld]) except +
        Integer num_steps()
        void step() except +
        bool step(Real) except +
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
        shared_ptr[Cpp_MesoscopicWorld] world()
        void run(Real) except +
        void run(Real, shared_ptr[Cpp_Observer]) except +
        void run(Real, vector[shared_ptr[Cpp_Observer]]) except +

## MesoscopicSimulator
#  a python wrapper for Cpp_MesoscopicSimulator
cdef class MesoscopicSimulator:
    cdef Cpp_MesoscopicSimulator* thisptr

cdef MesoscopicSimulator MesoscopicSimulator_from_Cpp_MesoscopicSimulator(Cpp_MesoscopicSimulator* s)

## Cpp_MesoscopicFactory
#  ecell4::meso::MesoscopicFactory
cdef extern from "ecell4/meso/MesoscopicFactory.hpp" namespace "ecell4::meso":
    cdef cppclass Cpp_MesoscopicFactory "ecell4::meso::MesoscopicFactory":
        Cpp_MesoscopicFactory(Cpp_Integer3&, Real) except +
        Cpp_MesoscopicWorld* create_world()
        Cpp_MesoscopicWorld* create_world(string)
        Cpp_MesoscopicWorld* create_world(Cpp_Real3&)
        Cpp_MesoscopicWorld* create_world(shared_ptr[Cpp_Model])
        Cpp_MesoscopicSimulator* create_simulator(shared_ptr[Cpp_Model], shared_ptr[Cpp_MesoscopicWorld])
        Cpp_MesoscopicSimulator* create_simulator(shared_ptr[Cpp_MesoscopicWorld])
        Cpp_MesoscopicFactory* rng_ptr(shared_ptr[Cpp_RandomNumberGenerator]&)
        @staticmethod
        Cpp_Integer3 default_matrix_sizes()
        @staticmethod
        Real default_subvolume_length()

## MesoscopicFactory
#  a python wrapper for Cpp_MesoscopicFactory
cdef class MesoscopicFactory:
    cdef Cpp_MesoscopicFactory* thisptr

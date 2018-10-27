from libcpp.string cimport string
from libcpp cimport bool

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *

from cpython cimport PyObject


cdef extern from "ecell4/ode/ODESimulator.hpp" namespace "ecell4::ode":
    cdef enum Cpp_ODESolverType "ecell4::ode::ODESolverType":
        Cpp_UNDEF "ecell4::ode::UNDEF"
        Cpp_RUNGE_KUTA_CASH_KARP54 "ecell4::ode::RUNGE_KUTA_CASH_KARP54"
        Cpp_ROSENBROCK4_CONTROLLER "ecell4::ode::ROSENBROCK4_CONTROLLER"
        Cpp_EULER "ecell4::ode::EULER"

## Cpp_ODEWorld
#  ecell4::ode::ODEWorld
cdef extern from "ecell4/ode/ODEWorld.hpp" namespace "ecell4::ode":
    cdef cppclass Cpp_ODEWorld "ecell4::ode::ODEWorld":
        Cpp_ODEWorld() except +
        Cpp_ODEWorld(Cpp_Real3&) except +
        Cpp_ODEWorld(string&) except +
        # SpaceTraits
        Real& t()
        void set_t(Real&)
        void reset(Cpp_Real3&)
        Cpp_Real3& edge_lengths()
        # CompartmentSpaceTraits
        Real &volume()
        Integer num_molecules(Cpp_Species &)
        Integer num_molecules_exact(Cpp_Species &)
        vector[Cpp_Species] list_species()

        # CompartmentSpace member functions
        void set_volume(Real &)
        pair[pair[Cpp_ParticleID, Cpp_Particle], bool] new_particle(Cpp_Particle& p)
        pair[pair[Cpp_ParticleID, Cpp_Particle], bool] new_particle(Cpp_Species& sp, Cpp_Real3& pos)
        void add_molecules(Cpp_Species &sp, Integer &num)
        void add_molecules(Cpp_Species &sp, Integer &num, shared_ptr[Cpp_Shape])
        void remove_molecules(Cpp_Species &sp, Integer &num)
        # Optional members
        Real get_value(Cpp_Species &)
        Real get_value_exact(Cpp_Species &)
        void set_value(Cpp_Species &sp, Real &num)
        void save(string) except +
        void load(string) except +
        bool has_species(Cpp_Species &)
        void reserve_species(Cpp_Species &)
        void release_species(Cpp_Species &)
        void bind_to(shared_ptr[Cpp_Model]) except +
        Real evaluate(Cpp_ReactionRule &) except +

## ODEWorld
#  a python wrapper for Cpp_ODEWorld
cdef class ODEWorld:
    cdef shared_ptr[Cpp_ODEWorld]* thisptr

cdef ODEWorld ODEWorld_from_Cpp_ODEWorld(shared_ptr[Cpp_ODEWorld] m)

## Cpp_ODESimulator
cdef extern from "ecell4/ode/ODESimulator.hpp" namespace "ecell4::ode":
    cdef cppclass Cpp_ODESimulator "ecell4::ode::ODESimulator":
        Cpp_ODESimulator(shared_ptr[Cpp_ODEWorld], shared_ptr[Cpp_Model], Cpp_ODESolverType) except+
        Cpp_ODESimulator(shared_ptr[Cpp_ODEWorld], shared_ptr[Cpp_Model]) except+
        Cpp_ODESimulator(shared_ptr[Cpp_ODEWorld], Cpp_ODESolverType) except+
        Cpp_ODESimulator(shared_ptr[Cpp_ODEWorld]) except+

        void initialize()
        void step() except +
        bool step(Real) except +
        Real next_time()
        Real t()
        void set_t(Real)
        Real dt()
        void set_dt(Real)
        Integer num_steps()
        bool check_reaction()
        Real absolute_tolerance() const
        Real relative_tolerance() const
        void set_absolute_tolerance(Real)
        void set_relative_tolerance(Real)

        shared_ptr[Cpp_Model] model()
        shared_ptr[Cpp_ODEWorld] world()

        void run(Real) except +
        void run(Real, shared_ptr[Cpp_Observer]) except +
        void run(Real, vector[shared_ptr[Cpp_Observer]]) except +
        void run(Real, is_dirty) except +
        void run(Real, shared_ptr[Cpp_Observer], is_dirty) except +
        void run(Real, vector[shared_ptr[Cpp_Observer]], is_dirty) except +

cdef class ODESimulator:
    cdef Cpp_ODESimulator *thisptr

cdef ODESimulator ODESimulator_from_Cpp_ODESimulator(Cpp_ODESimulator* s)

## Cpp_ODEFactory
#  ecell4::ode::ODEFactory
cdef extern from "ecell4/ode/ODEFactory.hpp" namespace "ecell4::ode":
    cdef cppclass Cpp_ODEFactory "ecell4::ode::ODEFactory":
        Cpp_ODEFactory(Cpp_ODESolverType, Real, Real, Real) except +
        Cpp_ODEFactory() except +
        Cpp_ODEWorld* world()
        Cpp_ODEWorld* world(string)
        Cpp_ODEWorld* world(Cpp_Real3&)
        Cpp_ODESimulator* simulator(shared_ptr[Cpp_ODEWorld], shared_ptr[Cpp_Model])
        Cpp_ODESimulator* simulator(shared_ptr[Cpp_ODEWorld])
        Cpp_ODEFactory* rng_ptr(shared_ptr[Cpp_RandomNumberGenerator]&)
        @staticmethod
        Cpp_ODESolverType default_solver_type()
        @staticmethod
        Real default_dt()
        @staticmethod
        Real default_abs_tol()
        @staticmethod
        Real default_rel_tol()

## ODEFactory
#  a python wrapper for Cpp_ODEFactory
cdef class ODEFactory:
    cdef Cpp_ODEFactory* thisptr

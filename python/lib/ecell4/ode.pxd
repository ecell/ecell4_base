from libcpp.string cimport string
from libcpp cimport bool

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *


## Cpp_ODEWorld
#  ecell4::ode::ODEWorld
cdef extern from "ecell4/ode/ODEWorld.hpp" namespace "ecell4::ode":
    cdef cppclass Cpp_ODEWorld "ecell4::ode::ODEWorld":
        Cpp_ODEWorld() except +
        Cpp_ODEWorld(Cpp_Position3&) except +
        Cpp_ODEWorld(string&) except +
        # SpaceTraits
        Real& t()
        void set_t(Real&)
        void reset(Cpp_Position3&)
        Cpp_Position3 edge_lengths()
        # CompartmentSpaceTraits
        Real &volume()
        Integer num_molecules(Cpp_Species &)
        Integer num_molecules_exact(Cpp_Species &)
        vector[Cpp_Species] list_species()

        # CompartmentSpace member functions
        void set_volume(Real &)
        void add_molecules(Cpp_Species &sp, Integer &num)
        void add_molecules(Cpp_Species &sp, Integer &num, Cpp_Shape&)
        void remove_molecules(Cpp_Species &sp, Integer &num)
        # Optional members
        Real get_value(Cpp_Species &)
        void set_value(Cpp_Species &sp, Real &num)
        void save(string)
        void load(string)
        bool has_species(Cpp_Species &)
        void reserve_species(Cpp_Species &)
        void release_species(Cpp_Species &)
        void bind_to(shared_ptr[Cpp_NetworkModel])

## ODEWorld
#  a python wrapper for Cpp_ODEWorld
cdef class ODEWorld:
    cdef shared_ptr[Cpp_ODEWorld]* thisptr

cdef ODEWorld ODEWorld_from_Cpp_ODEWorld(shared_ptr[Cpp_ODEWorld] m)

## Cpp_ODESimulator
#  ecell4::ode::ODESimulator
cdef extern from "ecell4/ode/ODESimulator.hpp" namespace "ecell4::ode":
    cdef cppclass Cpp_ODESimulator "ecell4::ode::ODESimulator":
        Cpp_ODESimulator(
            shared_ptr[Cpp_NetworkModel], shared_ptr[Cpp_ODEWorld]) except +
        void initialize()
        Real t()
        Integer num_steps()
        Real dt()
        Real next_time()
        void step()
        bool step(Real&)
        # Optional members
        void set_t(Real&)
        void set_dt(Real &)
        vector[Cpp_ReactionRule] last_reactions()
        shared_ptr[Cpp_NetworkModel] model()
        shared_ptr[Cpp_ODEWorld] world()
        void run(Real)
        void run(Real, vector[shared_ptr[Cpp_Observer]])

## ODESimulator
#  a python wrapper for Cpp_ODESimulator
cdef class ODESimulator:
    cdef Cpp_ODESimulator *thisptr

cdef ODESimulator ODESimulator_from_Cpp_ODESimulator(Cpp_ODESimulator* s)

## Cpp_ODEFactory
#  ecell4::ode::ODEFactory
cdef extern from "ecell4/ode/ODEFactory.hpp" namespace "ecell4::ode":
    cdef cppclass Cpp_ODEFactory "ecell4::ode::ODEFactory":
        Cpp_ODEFactory() except +
        Cpp_ODEWorld* create_world(string)
        Cpp_ODEWorld* create_world(Cpp_Position3&)
        Cpp_ODESimulator* create_simulator(shared_ptr[Cpp_Model], shared_ptr[Cpp_ODEWorld])
        Cpp_ODESimulator* create_simulator(shared_ptr[Cpp_ODEWorld])

## ODEFactory
#  a python wrapper for Cpp_ODEFactory
cdef class ODEFactory:
    cdef Cpp_ODEFactory* thisptr

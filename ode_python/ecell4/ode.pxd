from libcpp.string cimport string
from libcpp cimport bool

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *


## Cpp_ODEWorld
#  ecell4::ode::ODEWorld
cdef extern from "ecell4/ode/ODEWorld.hpp" namespace "ecell4::ode":
    cdef cppclass Cpp_ODEWorld "ecell4::ode::ODEWorld":
        Cpp_ODEWorld(Cpp_Position3&) except +
        # SpaceTraits
        Real& t()
        void set_t(Real&)
        void set_edge_lengths(Cpp_Position3&)
        Cpp_Position3 edge_lengths()
        # CompartmentSpaceTraits
        Real &volume()
        Real num_molecules(Cpp_Species &)
        vector[Cpp_Species] list_species()

        # CompartmentSpace member functions
        void set_volume(Real &)
        void add_molecules(Cpp_Species &sp, Real &num)
        void remove_molecules(Cpp_Species &sp, Real &num)
        # Optional members
        void set_num_molecules(Cpp_Species &sp, Real &num)
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

## ODESimulator
#  a python wrapper for Cpp_ODESimulator
cdef class ODESimulator:
    cdef Cpp_ODESimulator *thisptr

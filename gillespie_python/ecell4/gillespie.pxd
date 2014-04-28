from libcpp.string cimport string
from libcpp cimport bool

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *


## Cpp_GillespieWorld
#  ecell4::gillespie::GillespieWorld
cdef extern from "ecell4/gillespie/GillespieWorld.hpp" namespace "ecell4::gillespie":
    cdef cppclass Cpp_GillespieWorld "ecell4::gillespie::GillespieWorld":
        Cpp_GillespieWorld(Cpp_Position3&) except +
        Cpp_GillespieWorld(Cpp_Position3&, shared_ptr[Cpp_GSLRandomNumberGenerator]) except +
        void set_t(Real)
        Real t()
        Real volume()
        void set_edge_lengths(Cpp_Position3&)
        Cpp_Position3 edge_lengths()
        Integer num_molecules(Cpp_Species &)
        void add_molecules(Cpp_Species &sp, Integer &num)
        void remove_molecules(Cpp_Species &sp, Integer &num)
        void save(string)
        void load(string)
        void bind_to(shared_ptr[Cpp_NetworkModel])

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
            shared_ptr[Cpp_NetworkModel], shared_ptr[Cpp_GillespieWorld]) except +
        Integer num_steps()
        void step()
        bool step(Real)
        Real t()
        void set_t(Real)
        void set_dt(Real)
        Real dt()
        Real next_time()
        void initialize()
        Cpp_GSLRandomNumberGenerator& rng()
        shared_ptr[Cpp_NetworkModel] model()
        shared_ptr[Cpp_GillespieWorld] world()

## GillespieSimulator
#  a python wrapper for Cpp_GillespieSimulator
cdef class GillespieSimulator:
    cdef Cpp_GillespieSimulator* thisptr

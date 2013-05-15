from cython.operator cimport dereference as deref
from libcpp.string cimport string
from libcpp cimport bool

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *


## Cpp_GillespieWorld
#  ecell4::gillespie::GillespieWorld
cdef extern from "ecell4/gillespie/GillespieWorld.hpp" namespace "ecell4::gillespie":
    cdef cppclass Cpp_GillespieWorld "ecell4::gillespie::GillespieWorld":
        Cpp_GillespieWorld(Real, shared_ptr[Cpp_GSLRandomNumberGenerator]) except +
        void set_t(Real)
        Real t()
        Real volume()
        Integer num_species()
        bool has_species(Cpp_Species &)
        Integer num_molecules(Cpp_Species &)
        void add_species(Cpp_Species &)
        void remove_species(Cpp_Species &)
        void add_molecules(Cpp_Species &sp, Integer &num)
        void remove_molecules(Cpp_Species &sp, Integer &num)
        void save(string)

## GillespieWorld
#  a python wrapper for Cpp_GillespieWorld
cdef class GillespieWorld:
    cdef shared_ptr[Cpp_GillespieWorld]* thisptr

    def __cinit__(self, Real vol, GSLRandomNumberGenerator rng):
        # XXX: GSLRandomNumberGenerator -> RandomNumberGenerator
        self.thisptr = new shared_ptr[Cpp_GillespieWorld](
            new Cpp_GillespieWorld(vol, deref(rng.thisptr)))

    def __dealloc__(self):
        # XXX: Here, we release shared pointer,
        #      and if reference count to the GillespieWorld object,
        #      it will be released automatically.
        del self.thisptr

    def set_t(self, Real t):
        self.thisptr.get().set_t(t)

    def t(self):
        return self.thisptr.get().t()

    def volume(self):
        return self.thisptr.get().volume()

    def num_species(self):
        return self.thisptr.get().num_species()

    def has_species(self, Species sp):
        return self.thisptr.get().has_species(deref(sp.thisptr))

    def num_molecules(self, Species sp):
        return self.thisptr.get().num_molecules(deref(sp.thisptr))

    def add_species(self, Species sp):
        self.thisptr.get().add_species(deref(sp.thisptr))

    def remove_species(self, Species sp):
        self.thisptr.get().remove_species(deref(sp.thisptr))

    def add_molecules(self, Species sp, Integer num):
        self.thisptr.get().add_molecules(deref(sp.thisptr), num)

    def remove_species(self, Species sp, Integer num):
        self.thisptr.get().remove_molecules(deref(sp.thisptr), num)

    def save(self, string filename):
        self.thisptr.get().save(filename)

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
        Real dt()
        void initialize()
        Cpp_GSLRandomNumberGenerator& rng()

## GillespieSimulator
#  a python wrapper for Cpp_GillespieSimulator
cdef class GillespieSimulator:
    cdef Cpp_GillespieSimulator* thisptr

    def __cinit__(self, NetworkModel m, GillespieWorld w):
        self.thisptr = new Cpp_GillespieSimulator(
            deref(m.thisptr), deref(w.thisptr))

    def __dealloc__(self):
        del self.thisptr

    def num_steps(self):
        return self.thisptr.num_steps()

    def step(self):
        self.thisptr.step()

    def step_upto(self, Real upto):
        return self.thisptr.step(upto)

    def t(self):
        return self.thisptr.t()

    def set_t(self, Real new_t):
        self.thisptr.set_t(new_t)

    def initialize(self):
        self.thisptr.initialize()

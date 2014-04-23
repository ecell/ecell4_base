from cython.operator cimport dereference as deref
from libcpp.string cimport string

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *


## GillespieWorld
#  a python wrapper for Cpp_GillespieWorld
cdef class GillespieWorld:

    def __cinit__(self, Position3 edge_lengths, GSLRandomNumberGenerator rng = None):
        if rng is None:
            self.thisptr = new shared_ptr[Cpp_GillespieWorld](
                new Cpp_GillespieWorld(deref(edge_lengths.thisptr)))
        else:
            # XXX: GSLRandomNumberGenerator -> RandomNumberGenerator
            self.thisptr = new shared_ptr[Cpp_GillespieWorld](
                new Cpp_GillespieWorld(
                    deref(edge_lengths.thisptr), deref(rng.thisptr)))

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

    def num_molecules(self, Species sp):
        return self.thisptr.get().num_molecules(deref(sp.thisptr))

    def add_molecules(self, Species sp, Integer num):
        self.thisptr.get().add_molecules(deref(sp.thisptr), num)

    def remove_molecules(self, Species sp, Integer num):
        self.thisptr.get().remove_molecules(deref(sp.thisptr), num)

    def save(self, string filename):
        self.thisptr.get().save(filename)

    def load(self, string filename):
        self.thisptr.get().load(filename)

    def bind_to(self, NetworkModel m):
        self.thisptr.get().bind_to(deref(m.thisptr))

## GillespieSimulator
#  a python wrapper for Cpp_GillespieSimulator
cdef class GillespieSimulator:

    def __cinit__(self, NetworkModel m, GillespieWorld w):
        self.thisptr = new Cpp_GillespieSimulator(
            deref(m.thisptr), deref(w.thisptr))

    def __dealloc__(self):
        del self.thisptr

    def num_steps(self):
        return self.thisptr.num_steps()

    def step(self, upto = None):
        if upto is None:
            self.thisptr.step()
        else:
            return self.thisptr.step(<Real> upto)

    def t(self):
        return self.thisptr.t()

    def set_t(self, Real new_t):
        self.thisptr.set_t(new_t)

    def initialize(self):
        self.thisptr.initialize()

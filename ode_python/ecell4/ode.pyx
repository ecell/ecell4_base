from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address
from libcpp.string cimport string
from libcpp.vector cimport vector

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *

## ODEWorld
#  a python wrapper for Cpp_ODEWorld
cdef class ODEWorld:

    def __cinit__(self, Position3 edge_lengths):
        # XXX: GSLRandomNumberGenerator -> RandomNumberGenerator
        self.thisptr = new shared_ptr[Cpp_ODEWorld](
            new Cpp_ODEWorld(deref(edge_lengths.thisptr)) )

    def __dealloc__(self):
        # XXX: Here, we release shared pointer,
        #      and if reference count to the ODEWorld object become zero,
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

    def list_species(self):
        cdef vector[Cpp_Species] raw_list_species = self.thisptr.get().list_species() 
        retval = []
        cdef vector[Cpp_Species].iterator it = raw_list_species.begin()
        while it != raw_list_species.end():
            retval.append( 
                Species_from_Cpp_Species(<Cpp_Species*> (address(deref(it))))) 
            inc(it)
        return retval

    def set_volume(self, Real vol):
        self.thisptr.get().set_volume(vol)

    def add_molecules(self, Species sp, Real num):
        self.thisptr.get().add_molecules(deref(sp.thisptr), num)

    def remove_species(self, Species sp, Real num):
        self.thisptr.get().remove_molecules(deref(sp.thisptr), num)

    def set_num_molecules(self, Species sp, Real num):
        self.thisptr.get().set_num_molecules(deref(sp.thisptr), num)

    def save(self, string filename):
        self.thisptr.get().save(filename)

    def has_species(self, Species sp):
        return self.thisptr.get().has_species(deref(sp.thisptr))

    def reserve_species(self, Species sp):
        self.thisptr.get().reserve_species(deref(sp.thisptr))

    def release_species(self, Species sp):
        self.thisptr.get().release_species(deref(sp.thisptr))

## ODESimulator
#  a python wrapper for Cpp_ODESimulator
cdef class ODESimulator:

    def __cinit__(self, NetworkModel m, ODEWorld w):
        self.thisptr = new Cpp_ODESimulator(
            deref(m.thisptr), deref(w.thisptr))

    def __dealloc__(self):
        del self.thisptr

    def num_steps(self):
        return self.thisptr.num_steps()

    def step(self, upto = None):
        if upto is None:
            self.thisptr.step()
        else:
            return self.thisptr.step(upto)

    def t(self):
        return self.thisptr.t()

    def set_t(self, Real new_t):
        self.thisptr.set_t(new_t)

    def initialize(self):
        self.thisptr.initialize()



from cython.operator cimport dereference as deref

from libcpp.vector cimport vector
from libcpp cimport set
from libcpp cimport bool

include "types.pxi"

from PyEcell4 cimport *
from PyEcell4 import PySpecies


cdef extern from "ecell4/gillespie/GillespieWorld.hpp" namespace "ecell4::gillespie":
    cdef cppclass GillespieWorld:
        GillespieWorld(Real) except +
        void set_t(Real)
        Real t()
        Real volume()
        Integer num_species()
        bool has_species(Species &)
        Integer num_molecules(Species &)
        
        void add_species(Species &)
        void remove_species(Species &)
        void add_molecules(Species &sp, Integer &num)
        void remove_molecules(Species &sp, Integer &num)

cdef class PyGillespieWorld:
    cdef GillespieWorld *thisptr
    def __cinit__(self, Real vol):
        self.thisptr = new GillespieWorld(vol)
    def __dealloc__(self):
        del self.thisptr
    
    def set_t(self, Real t):
        self.thisptr.set_t(t)
    def t(self):
        return self.thisptr.t()
    def volume(self):
        return self.thisptr.volume()
    def num_species(self):
        return self.thisptr.num_species()
    def has_species(self, PySpecies sp):
        return self.thisptr.has_species( deref(sp.thisptr) )
    def num_molecules(self, PySpecies sp):
        return self.thisptr.num_molecules( deref(sp.thisptr) )

    def add_species(self, PySpecies sp):
        self.thisptr.add_species(deref(sp.thisptr) )
    def remove_species(self, PySpecies sp):
        self.thisptr.remove_species(deref(sp.thisptr))
    def add_molecules(self, PySpecies sp, Integer num):
        self.thisptr.add_molecules(deref(sp.thisptr), num)
    def remove_species(self, PySpecies sp, Integer num):
        self.thisptr.remove_molecules(deref(sp.thisptr), num)


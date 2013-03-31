# distutils: language = c++
# distutils: sources = ../core/CompartmentSpace.cpp

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string 
from cython.operator cimport dereference as deref

from PySpecies cimport Species, PySpecies
from PySpecies import PySpecies

#ctypedef double Real 
#ctypedef int Integer

cdef extern from "ecell4/core/CompartmentSpace.hpp" namespace "ecell4":
    cdef cppclass CompartmentSpaceVectorImpl:
        #Constructor
        CompartmentSpaceVectorImpl(Real) except+
        Real volume()
        Integer num_species()
        bool has_species(Species &sp)
        Integer num_molecules(Species &sp)
        void set_volume(Real)
        void add_species(Species &sp)
        void remove_species(Species &sp)
        void add_molecules(Species &sp, Integer num)
        void remove_molecules(Species &sp, Integer num)

cdef class PyCompartmentSpace:
    #cdef CompartmentSpaceVectorImpl *thisptr
    def __cinit__(self, Real volume):
        self.thisptr = new CompartmentSpaceVectorImpl(volume)
    def __dealloc__(self):
        del self.thisptr
    def volume(self):
        return self.thisptr.volume()

    def num_species(self):
        return self.thisptr.num_species()

    def has_species(self, PySpecies sp):
        return self.thisptr.has_species( deref(sp.thisptr) )

    def num_molecules(self, PySpecies sp):
        return self.thisptr.num_molecules( deref(sp.thisptr) )
    def set_volume(self, Real volume):
        self.thisptr.set_volume(volume)
    def add_species(self, PySpecies sp):
        self.thisptr.add_species( deref(sp.thisptr) )
    def remove_species(self, PySpecies sp):
        self.thisptr.remove_species( deref(sp.thisptr) )
    def add_molecules(self, PySpecies sp, Integer num):
        self.thisptr.add_molecules( deref(sp.thisptr), num )
    def remove_molecules(self, PySpecies sp, Integer num):
        self.thisptr.remove_molecules( deref(sp.thisptr), num )

# distutils: language = c++
# distutils: sources = ../ode/ODEWorld.cpp

from libcpp cimport bool
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref

ctypedef double Real

#cimport PySpecies;  
from PySpecies cimport Species,PySpecies
from PySpecies import PySpecies

cdef extern from "ecell4/ode/ODEWorld.hpp" namespace "ecell4::ode":
    cdef cppclass ODEWorld:
        ODEWorld(Real volume) except +
        Real & volume()
        void set_volume(Real & volume)
        int num_species()
        bool has_species(Species sp)
        void add_species(Species &sp) 
        void remove_species(Species sp)
        Real num_molecules(Species sp)
        void set_num_molecules(Species sp, Real num)
        void add_molecules(Species sp, Real num) 
        void remove_molecules(Species sp, Real num)
        
        
cdef class PyOdeWorld:
    cdef ODEWorld *thisptr
    def __cinit__(self, Real vol):
        self.thisptr = new ODEWorld(vol)
    def __dealloc__(self):
        del self.thisptr

    def volume(self):
        return self.thisptr.volume()

    def set_volume(self, Real vol):
        self.thisptr.set_volume(vol)

    def num_species(self):
        return self.thisptr.num_species()

    def has_species(self, PySpecies sp):
        return self.thisptr.has_species( deref(sp.thisptr) )

    def add_species(self, PySpecies sp):
        self.thisptr.add_species( deref(sp.thisptr) ) 

    def remove_species(self, PySpecies sp):
        self.thisptr.remove_species( deref(sp.thisptr) )
    def num_molecules(self, PySpecies sp):
        return self.thisptr.num_molecules( deref(sp.thisptr) )
    def set_num_molecules(self, PySpecies sp, Real num):
        self.thisptr.set_num_molecules( deref(sp.thisptr), num)
    def add_molecules(self, PySpecies sp, Real num):
        self.thisptr.add_molecules( deref(sp.thisptr), num )
    def remove_molecules(self, PySpecies sp, Real num):
        self.thisptr.remove_molecules( deref(sp.thisptr), num )



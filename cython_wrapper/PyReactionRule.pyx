# distutils: language = c++
# distutils: sources = ../core/ReactionRule.cpp

from cython.operator cimport dereference as deref

from libcpp.vector cimport vector
from libcpp cimport set
from libcpp cimport bool


from PySpecies cimport Species,PySpecies
from PySpecies import PySpecies

cdef extern from "ecell4/core/ReactionRule.hpp" namespace "ecell4":
    cdef cppclass ReactionRule:
        ReactionRule() except +
        Real k()
        #set[Species] reactants() 
        #set[Species] products()
        void set_k(Real)
        void add_reactant(Species)
        void add_product(Species)

cdef class PyReactionRule:
 #   cdef ReactionRule *thisptr
    def __cinit__(self):
        self.thisptr = new ReactionRule() 
    def __dealloc__(self):
        del self.thisptr

    def k(self):
        return self.thisptr.k()

    def set_k(self, Real k):
        self.thisptr.set_k(k)

    def add_reactant(self, PySpecies sp):
        self.thisptr.add_reactant( deref(sp.thisptr) )

    def add_product(self, PySpecies sp):
        self.thisptr.add_product( deref(sp.thisptr) )


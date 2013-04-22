
from cython.operator cimport dereference as deref
from libcpp.vector cimport vector
from libcpp cimport set
from libcpp cimport bool


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

def create_unimolecular_reaction_rule(
        PySpecies reactant1, PySpecies product1, Real k):
    rr = PyReactionRule()
    rr.add_reactant(reactant1)
    rr.add_product(product1)
    rr.set_k(k)
    return rr

def create_binding_reaction_rule(
        PySpecies reactant1, PySpecies reactant2, PySpecies product1, Real k):
    rr = PyReactionRule()
    rr.add_reactant(reactant1)
    rr.add_reactant(reactant2)
    rr.add_product(product1)
    rr.set_k(k)
    return rr

def create_unbinding_reaction_rue(
        PySpecies reactant1, PySpecies product1, PySpecies product2, Real k):
    rr = PyReactionRule()
    rr.add_reactant(reactant1)
    rr.add_product(product1)
    rr.add_product(product2)
    rr.set_k(k)
    return rr

def create_repulsive_reaction_rule(
        PySpecies reactant1, PySpecies reactant2):
    rr = PyReactionRule()
    rr.set_k(0.0)
    rr.add_reactant(reactant1)
    rr.add_product(reactant2)
    return rr

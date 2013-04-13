# distutils: language = c++
# distutils: sources = ../core/NetworkModel.cpp

from cython.operator cimport dereference as deref

from libcpp.vector cimport vector
from libcpp cimport set
from libcpp cimport bool

include "types.pxi"

from PySpecies cimport Species,PySpecies
from PySpecies import PySpecies

from PyReactionRule cimport ReactionRule, PyReactionRule
from PyReactionRule import PyReactionRule

cdef extern from "ecell4/core/NetworkModel.hpp" namespace "ecell4":
    cdef cppclass NetworkModel:
        NetworkModel() except +
        void add_species(Species sp)
        bool has_species(Species sp)
        void remove_species(Species sp)
        void add_reaction_rule(ReactionRule)
        void remove_reaction_rule(ReactionRule)
        bool has_reaction_rule(ReactionRule)

cdef class PyNetworkModel:
    cdef NetworkModel *thisptr
    def __cinit__(self):
        self.thisptr = new NetworkModel() 
    def __dealloc__(self):
        del self.thisptr
    #def k(self):
    #    return self.thisptr.k()
    
    # HANDLERS FOR SPECIES
    def add_species(self, PySpecies sp):
        self.thisptr.add_species( deref(sp.thisptr) )
    def has_species(self, PySpecies sp):
        return self.thisptr.has_species( deref(sp.thisptr) )
    def remove_species(self, PySpecies sp):
        self.thisptr.remove_species( deref(sp.thisptr) )

    # HANDLERS FOR REACTION_RULES
    def add_reaction_rule(self, PyReactionRule rr):
        self.thisptr.add_reaction_rule( deref(rr.thisptr) )
    def remove_reaction_rule(self, PyReactionRule rr):
        self.thisptr.remove_reaction_rule( deref(rr.thisptr) )
    def has_reaction_rule(self, PyReactionRule rr):
        self.thisptr.has_reaction_rule( deref(rr.thisptr) )
    
'''
    def add_reactant(self, PySpecies sp):
        self.thisptr.add_reactant( deref(sp.thisptr) )
    def add_product(self, PySpecies sp):
        self.thisptr.add_product( deref(sp.thisptr) )
    def reactants(self):
        #self.thisptr.reactants() 
        pass
'''

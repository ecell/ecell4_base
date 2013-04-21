
from cython.operator cimport dereference as deref

from libcpp.vector cimport vector
from libcpp cimport set
from libcpp cimport bool

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
    #cdef NetworkModel *thisptr
    def __cinit__(self):
        #self.thisptr = new NetworkModel() 
        self.thisptr = new shared_ptr[NetworkModel](new NetworkModel())
    def __dealloc__(self):
        del self.thisptr
    # HANDLERS FOR SPECIES
    def add_species(self, PySpecies sp):
        self.thisptr.get().add_species( deref(sp.thisptr) )
    def has_species(self, PySpecies sp):
        return self.thisptr.get().has_species( deref(sp.thisptr) )
    def remove_species(self, PySpecies sp):
        self.thisptr.get().remove_species( deref(sp.thisptr) )

    # HANDLERS FOR REACTION_RULES
    def add_reaction_rule(self, PyReactionRule rr):
        self.thisptr.get().add_reaction_rule( deref(rr.thisptr) )
    def remove_reaction_rule(self, PyReactionRule rr):
        self.thisptr.get().remove_reaction_rule( deref(rr.thisptr) )
    def has_reaction_rule(self, PyReactionRule rr):
        self.thisptr.get().has_reaction_rule( deref(rr.thisptr) )
    
'''
    def add_reactant(self, PySpecies sp):
        self.thisptr.add_reactant( deref(sp.thisptr) )
    def add_product(self, PySpecies sp):
        self.thisptr.add_product( deref(sp.thisptr) )
    def reactants(self):
        #self.thisptr.reactants() 
        pass
'''

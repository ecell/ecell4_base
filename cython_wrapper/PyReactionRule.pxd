

from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp cimport bool

include "types.pxi"

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
    cdef ReactionRule *thisptr


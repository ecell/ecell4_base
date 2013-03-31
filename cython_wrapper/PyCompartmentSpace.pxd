
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string 
from cython.operator cimport dereference as deref

from PySpecies cimport Species, PySpecies
from PySpecies import PySpecies

ctypedef double Real 
ctypedef int Integer

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
    cdef CompartmentSpaceVectorImpl *thisptr


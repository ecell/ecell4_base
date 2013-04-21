
#============================================================
#   Common Declaration and imports
#============================================================

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.set cimport set
from libcpp cimport bool

include "types.pxi"


#============================================================
#   Boost.shared_ptr<T>
#============================================================
cdef extern from "boost/shared_ptr.hpp" namespace "boost":
    cdef cppclass shared_ptr[T]:
        shared_ptr(T *ptr)
        T* get()

#============================================================
#   Species
#============================================================
cdef extern from "ecell4/core/Species.hpp" namespace "ecell4":
    cdef cppclass Species:
        # Constructor
        Species(string) except +
        Species(string, string)
        Species(string, string, string)
#       serial_type serial()
        string name()
        string get_attribute(string)
        void set_attribute(string,string)
        void remove_attribute(string)

cdef class PySpecies:
    cdef Species *thisptr

#============================================================
#   ReactionRule
#============================================================
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

#============================================================
#   CompartmentSpace
#============================================================
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

#============================================================
#   NetworkModel
#============================================================
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

#============================================================
#   Position3
#============================================================
cdef extern from "ecell4/core/Position3.hpp" namespace "ecell4":
    cdef cppclass Position3:
        Position3() except +
        Position3(Real, Real, Real) except +
        Position3(Position3 &rhs) except+
    Position3 add(Position3, Position3)
    Position3 subtract(Position3, Position3)
    Position3 divide(Position3, Real)
    Position3 multiply(Position3, Real)
    Position3 modulo(Position3, Real)
    Position3 modulo(Position3, Position3)
    Position3 abs(Position3)
    Real dot_product(Position3, Position3)
    Position3 cross_product(Position3, Position3)
    Integer length_sq(Position3)
    Integer length(Position3)
    Position3 operator+(Position3, Position3)
    Position3 operator-(Position3, Position3)
    Position3 operator/(Position3, Real)
    Position3 operator*(Position3, Real)

cdef class PyPosition3:
    cdef Position3 *thisptr



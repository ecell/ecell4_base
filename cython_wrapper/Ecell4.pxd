
#============================================================
#   Common Declaration and imports
#============================================================

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.set cimport set
from libcpp cimport bool

include "types.pxi"

#============================================================
#   Stl:    MultiSet
#============================================================
cdef extern from "<set>" namespace "std":
    cdef cppclass multiset[T]:
        multiset() except +
        multiset(multiset &) except+
        cppclass iterator:
            T& operator*()
            iterator operator++() 
            iterator operator--() 
            bint operator==(iterator) 
            bint operator!=(iterator) 
        iterator begin() 
        iterator end() 

#============================================================
#   Boost.shared_ptr<T>
#============================================================
cdef extern from "<boost/shared_ptr.hpp>" namespace "boost":
    cdef cppclass shared_ptr[T]:
        shared_ptr(T *ptr)
        T* get()

#============================================================
#   RandomNumberGenerator
#============================================================
cdef extern from "gsl/gsl_rng.h":
    ctypedef struct gsl_rng:
        pass

cdef extern from "ecell4/core/RandomNumberGenerator.hpp" namespace "ecell4":
    cdef cppclass Cpp_GSLRandomNumberGenerator "ecell4::GSLRandomNumberGenerator":
        #GSLRandomNumberGenerator(shared_ptr[gsl_rng]) except +
        Cpp_GSLRandomNumberGenerator() except +
        Real uniform(Real, Real)
        Integer uniform_int(Integer, Integer)
        Real gaussian(Real, Real)
        void seed(Integer)

# For Dereference
cdef class RandomNumberGenerator:
    cdef Cpp_GSLRandomNumberGenerator *thisptr

#============================================================
#   Species
#============================================================
cdef extern from "ecell4/core/Species.hpp" namespace "ecell4":
    cdef cppclass Cpp_Species "ecell4::Species":
        # Constructor
        Cpp_Species(string) except +
        Cpp_Species(string, string) except +
        Cpp_Species(string, string, string) except +
        Cpp_Species(Cpp_Species &) except+
#       serial_type serial()
        string name()
        string get_attribute(string)
        void set_attribute(string,string)
        void remove_attribute(string)

cdef class Species:
    cdef Cpp_Species *thisptr

#============================================================
#   ReactionRule
#============================================================
cdef extern from "ecell4/core/ReactionRule.hpp" namespace "ecell4":
    cdef cppclass Cpp_ReactionRule "ecell4::ReactionRule":
        Cpp_ReactionRule() except +
        Real k()
        multiset[Cpp_Species]& reactants() 
        multiset[Cpp_Species]& products()
        void set_k(Real)
        void add_reactant(Cpp_Species)
        void add_product(Cpp_Species)

cdef class ReactionRule:
    cdef Cpp_ReactionRule *thisptr

#============================================================
#   CompartmentSpace
#============================================================
cdef extern from "ecell4/core/CompartmentSpace.hpp" namespace "ecell4":
    cdef cppclass Cpp_CompartmentSpaceVector "ecell4::CompartmentSpaceVectorImpl":
        #Constructor
        Cpp_CompartmentSpaceVector(Real) except+
        Real volume()
        Integer num_species()
        bool has_species(Cpp_Species &sp)
        Integer num_molecules(Cpp_Species &sp)
        void set_volume(Real)
        void add_species(Cpp_Species &sp)
        void remove_species(Cpp_Species &sp)
        void add_molecules(Cpp_Species &sp, Integer num)
        void remove_molecules(Cpp_Species &sp, Integer num)

cdef class CompartmentSpace:
    cdef Cpp_CompartmentSpaceVector *thisptr

#============================================================
#   NetworkModel
#============================================================
cdef extern from "ecell4/core/NetworkModel.hpp" namespace "ecell4":
    cdef cppclass Cpp_NetworkModel "ecell4::NetworkModel":
        Cpp_NetworkModel() except +
        void add_species(Cpp_Species sp)
        bool has_species(Cpp_Species sp)
        void remove_species(Cpp_Species sp)
        void add_reaction_rule(Cpp_ReactionRule)
        void remove_reaction_rule(Cpp_ReactionRule)
        bool has_reaction_rule(Cpp_ReactionRule)

cdef class NetworkModel:
    #cdef Cpp_NetworkModel *thisptr
    cdef shared_ptr[Cpp_NetworkModel] *thisptr

#============================================================
#   Position3
#============================================================
cdef extern from "ecell4/core/Position3.hpp" namespace "ecell4":
    cdef cppclass Cpp_Position3 "ecell4::Position3":
        Cpp_Position3() except +
        Cpp_Position3(Real, Real, Real) except +
        Cpp_Position3(Cpp_Position3 &rhs) except+
    Cpp_Position3 add(Cpp_Position3, Cpp_Position3)
    Cpp_Position3 subtract(Cpp_Position3, Cpp_Position3)
    Cpp_Position3 divide(Cpp_Position3, Real)
    Cpp_Position3 multiply(Cpp_Position3, Real)
    Cpp_Position3 modulo(Cpp_Position3, Real)
    Cpp_Position3 modulo(Cpp_Position3, Cpp_Position3)
    Cpp_Position3 abs(Cpp_Position3)
    Real dot_product(Cpp_Position3, Cpp_Position3)
    Cpp_Position3 cross_product(Cpp_Position3, Cpp_Position3)
    Integer length_sq(Cpp_Position3)
    Integer length(Cpp_Position3)
    Cpp_Position3 operator+(Cpp_Position3, Cpp_Position3)
    Cpp_Position3 operator-(Cpp_Position3, Cpp_Position3)
    Cpp_Position3 operator/(Cpp_Position3, Real)
    Cpp_Position3 operator*(Cpp_Position3, Real)

cdef class Position3:
    cdef Cpp_Position3 *thisptr

#============================================================
#   Particle
#============================================================
cdef extern from "ecell4/core/Particle.hpp" namespace "ecell4":
    cdef cppclass Cpp_Particle "ecell4::Particle":
        Cpp_Particle() except +
        Cpp_Particle(Cpp_Species, Cpp_Position3, Real radius, Real D) except +
        Cpp_Position3 position()
        Real radius()
        Real D()
        Cpp_Species &species()

cdef class Particle:
    cdef Cpp_Particle *thisptr

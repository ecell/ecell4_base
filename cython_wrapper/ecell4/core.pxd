from libcpp.string cimport string
from libcpp cimport bool

# XXX Tomplorary using cython stl support.
#        Perhaps, we should consider importing std::pair by ourselves
#        that don't cast c-objects into python objects automatically.
from libcpp.pair cimport pair

from types cimport *
from multiset cimport multiset
from shared_ptr cimport shared_ptr


cdef extern from "gsl/gsl_rng.h":
    ctypedef struct gsl_rng:
        pass

## Cpp_GSLRandomNumberGenerator
#  ecell4::GSLRandomNumberGenerator
cdef extern from "ecell4/core/RandomNumberGenerator.hpp" namespace "ecell4":
    cdef cppclass Cpp_GSLRandomNumberGenerator "ecell4::GSLRandomNumberGenerator":
        # GSLRandomNumberGenerator(shared_ptr[gsl_rng]) except +
        Cpp_GSLRandomNumberGenerator() except +
        Real uniform(Real, Real)
        Integer uniform_int(Integer, Integer)
        Real gaussian(Real, Real)
        void seed(Integer)

## RandomNumberGenerator
#  a python wrapper for Cpp_GSLRandomNumberGenerator
cdef class GSLRandomNumberGenerator:
    # cdef Cpp_GSLRandomNumberGenerator* thisptr
    cdef shared_ptr[Cpp_GSLRandomNumberGenerator]* thisptr

## Cpp_Species
#  ecell4::Species
cdef extern from "ecell4/core/Species.hpp" namespace "ecell4":
    cdef cppclass Cpp_Species "ecell4::Species":
        Cpp_Species(string) except +
        Cpp_Species(string, string) except +
        Cpp_Species(string, string, string) except +
        Cpp_Species(Cpp_Species &) except+
        # serial_type serial()
        string name()
        string get_attribute(string)
        void set_attribute(string,string)
        void remove_attribute(string)

## Species
#  a python wrapper for Cpp_Species
cdef class Species:
    cdef Cpp_Species* thisptr

cdef Species Cpp_Species_to_Species(Cpp_Species *sp)

## Cpp_ReactionRule
#  ecell4::ReactionRule
cdef extern from "ecell4/core/ReactionRule.hpp" namespace "ecell4":
    cdef cppclass Cpp_ReactionRule "ecell4::ReactionRule":
        Cpp_ReactionRule() except +
        Real k()
        multiset[Cpp_Species]& reactants()
        multiset[Cpp_Species]& products()
        void set_k(Real)
        void add_reactant(Cpp_Species)
        void add_product(Cpp_Species)

## ReactionRule
#  a python wrapper for Cpp_ReactionRule
cdef class ReactionRule:
    cdef Cpp_ReactionRule* thisptr

## Cpp_CompartmentSpaceVectorImpl
#  ecell4::CompartmentSpaceVectorImpl
cdef extern from "ecell4/core/CompartmentSpace.hpp" namespace "ecell4":
    cdef cppclass Cpp_CompartmentSpaceVectorImpl "ecell4::CompartmentSpaceVectorImpl":
        Cpp_CompartmentSpaceVectorImpl(Real) except+
        Real volume()
        Integer num_species()
        bool has_species(Cpp_Species &sp)
        Integer num_molecules(Cpp_Species &sp)
        void set_volume(Real)
        void add_species(Cpp_Species &sp)
        void remove_species(Cpp_Species &sp)
        void add_molecules(Cpp_Species &sp, Integer num)
        void remove_molecules(Cpp_Species &sp, Integer num)

## CompartmentSpaceVectorImpl
#  a python wrapper for Cpp_CompartmentSpaceVectorImpl
cdef class CompartmentSpaceVectorImpl:
    cdef Cpp_CompartmentSpaceVectorImpl* thisptr

## Cpp_NetworkModel
#  ecell4::NetworkModel
cdef extern from "ecell4/core/NetworkModel.hpp" namespace "ecell4":
    cdef cppclass Cpp_NetworkModel "ecell4::NetworkModel":
        Cpp_NetworkModel() except +
        void add_species(Cpp_Species sp)
        bool has_species(Cpp_Species sp)
        void remove_species(Cpp_Species sp)
        void add_reaction_rule(Cpp_ReactionRule)
        void remove_reaction_rule(Cpp_ReactionRule)
        bool has_reaction_rule(Cpp_ReactionRule)

## NetworkModel
#  a python wrapper for Cpp_NetowrkModel, but wrapped by shared_ptr
cdef class NetworkModel:
    # cdef Cpp_NetworkModel* thisptr
    cdef shared_ptr[Cpp_NetworkModel]* thisptr

## Cpp_Position3
#  ecell4::Position3
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

## Position3
#  a python wrapper for Cpp_Position3
cdef class Position3:
    cdef Cpp_Position3* thisptr

cdef Position3 Cpp_Position3_to_Position3(Cpp_Position3 *p)

## Cpp_ParticleID
#  ecell4::ParticleID
cdef extern from "ecell4/core/Particle.hpp" namespace "ecell4":
    ctypedef Tbase_ "ecell4::ParticleID"
    ctypedef int lot_type
    ctypedef unsigned long long serial_type
    ctypedef pair[lot_type, serial_type] value_type

    cdef cppclass Cpp_ParticleID "ecell4::ParticleID":
        Cpp_ParticleID(value_type)
        Cpp_ParticleID log_add(lot_type &rhs)
        Cpp_ParticleID log_subtract(lot_type &rhs)
        Cpp_ParticleID &lot_advance(lot_type &rhs)
        Cpp_ParticleID &lot_retraace(lot_type &rhs)
        Cpp_ParticleID serial_add(serial_type &rhs)
        Cpp_ParticleID serial_subtract(serial_type &rhs)
        Cpp_ParticleID &serial_advance(serial_type &rhs)
        Cpp_ParticleID &serial_retrace(serial_type &rhs)
        # Cpp_ParticleID &operator=(Cpp_ParticleID &rhs) # XXX not yet suppoted
        bool operator==(Cpp_ParticleID &rhs)
        bool operator!=(Cpp_ParticleID &rhs)
        bool operator<(Cpp_ParticleID &rhs)
        bool operator>=(Cpp_ParticleID &rhs)
        bool operator>(Cpp_ParticleID &rhs)
        bool operator<=(Cpp_ParticleID &rhs)
        # operator value_type()
        value_type &operator() ()
        lot_type &lot()
        serial_type &serial()

## Cpp_Particle
#  ecell4::Particle
cdef extern from "ecell4/core/Particle.hpp" namespace "ecell4":
    cdef cppclass Cpp_Particle "ecell4::Particle":
        Cpp_Particle() except +
        Cpp_Particle(Cpp_Species, Cpp_Position3, Real radius, Real D) except +
        Cpp_Position3 position()
        Real radius()
        Real D()
        Cpp_Species &species()

## Particle
#  a python wrapper for Cpp_Particle
cdef class Particle:
    cdef Cpp_Particle* thisptr

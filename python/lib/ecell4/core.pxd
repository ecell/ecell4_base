from libcpp.string cimport string
from libcpp cimport bool

# XXX Tomplorary using cython stl support.
#        Perhaps, we should consider importing std::pair by ourselves
#        that don't cast c-objects into python objects automatically.
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp.map cimport map

from types cimport Real, Integer
from multiset cimport multiset
from shared_ptr cimport shared_ptr


cdef string tostring(ustr)

cdef extern from "gsl/gsl_rng.h":
    ctypedef struct gsl_rng:
        pass

## Cpp_GSLRandomNumberGenerator
#  ecell4::GSLRandomNumberGenerator
cdef extern from "ecell4/core/RandomNumberGenerator.hpp" namespace "ecell4":
    cdef cppclass Cpp_RandomNumberGenerator "ecell4::RandomNumberGenerator":
        # RandomNumberGenerator(shared_ptr[gsl_rng]) except +
        # Cpp_RandomNumberGenerator() except +
        Real random()
        Real uniform(Real, Real)
        Integer uniform_int(Integer, Integer)
        Real gaussian(Real, Real)
        Real gaussian(Real)
        Integer binomial(Real, Integer)
        void seed(Integer)
        void seed()
        void save(string) except +
        void load(string) except +

    cdef cppclass Cpp_GSLRandomNumberGenerator "ecell4::GSLRandomNumberGenerator":
        # GSLRandomNumberGenerator(shared_ptr[gsl_rng]) except +
        Cpp_GSLRandomNumberGenerator() except +
        Cpp_GSLRandomNumberGenerator(Integer) except +
        Cpp_GSLRandomNumberGenerator(string) except +
        Real uniform(Real, Real)
        Integer uniform_int(Integer, Integer)
        Real gaussian(Real, Real)
        Real gaussian(Real)
        void seed(Integer)
        void seed()
        void save(string) except +
        void load(string) except +

## RandomNumberGenerator
#  a python wrapper for Cpp_GSLRandomNumberGenerator
cdef class GSLRandomNumberGenerator:
    # cdef Cpp_GSLRandomNumberGenerator* thisptr
    # cdef shared_ptr[Cpp_GSLRandomNumberGenerator]* thisptr
    cdef shared_ptr[Cpp_RandomNumberGenerator]* thisptr

cdef GSLRandomNumberGenerator GSLRandomNumberGenerator_from_Cpp_RandomNumberGenerator(
    shared_ptr[Cpp_RandomNumberGenerator])

## Cpp_UnitSpecies
#  ecell4::UnitSpecies
cdef extern from "ecell4/core/UnitSpecies.hpp" namespace "ecell4":
    cdef cppclass Cpp_UnitSpecies "ecell4::UnitSpecies":
        Cpp_UnitSpecies() except +
        Cpp_UnitSpecies(string) except +
        Cpp_UnitSpecies(Cpp_UnitSpecies&) except+
        bool operator==(Cpp_UnitSpecies& rhs)
        bool operator<(Cpp_UnitSpecies& rhs)
        bool operator>(Cpp_UnitSpecies& rhs)
        string serial()
        string name()
        void deserialize(string) except+
        bool add_site(string, string, string)

## UnitSpecies
#  a python wrapper for Cpp_UnitSpecies
cdef class UnitSpecies:
    cdef Cpp_UnitSpecies* thisptr

cdef UnitSpecies UnitSpecies_from_Cpp_UnitSpecies(Cpp_UnitSpecies *sp)

cdef extern from "boost/variant.hpp" namespace "boost":
    cdef cppclass boost_variant "boost::variant" [T1, T2, T3, T4]:
        pass

    U* boost_get "boost::get" [U, T1, T2, T3, T4] (boost_variant[T1, T2, T3, T4]*) except +

ctypedef boost_variant[string, Real, Integer, bool] Cpp_Species_value_type

cdef boost_get_from_Cpp_Species_value_type(Cpp_Species_value_type value)

## Cpp_Species
#  ecell4::Species
cdef extern from "ecell4/core/Species.hpp" namespace "ecell4":
    cdef cppclass Cpp_Species "ecell4::Species":
        Cpp_Species() except +
        Cpp_Species(string) except +
        Cpp_Species(string, Real, Real) except +
        Cpp_Species(string, Real, Real, string) except +
        # Cpp_Species(string, string) except +
        Cpp_Species(string, string, string) except +
        Cpp_Species(string, string, string, string) except +
        Cpp_Species(Cpp_Species&) except+
        bool operator==(Cpp_Species& rhs)
        bool operator<(Cpp_Species& rhs)
        bool operator>(Cpp_Species& rhs)
        string serial() # string == serial_type
        # string get_attribute(string) except +
        Cpp_Species_value_type get_attribute(string) except +
        void set_attribute[T](string, T&)
        Integer count(Cpp_Species& sp) except +
        void remove_attribute(string) except +
        bool has_attribute(string)
        vector[pair[string, Cpp_Species_value_type]] list_attributes()
        void add_unit(Cpp_UnitSpecies)
        vector[Cpp_UnitSpecies]& units()
        Cpp_Species* D_ptr(string)
        Cpp_Species* radius_ptr(string)
        Cpp_Species* location_ptr(string)

## Species
#  a python wrapper for Cpp_Species
cdef class Species:
    cdef Cpp_Species* thisptr

cdef Species Species_from_Cpp_Species(Cpp_Species *sp)

## Cpp_ReactionRule
#  ecell4::ReactionRule
cdef extern from "ecell4/core/ReactionRule.hpp" namespace "ecell4":
    cdef enum Cpp_ReactionRulePolicyType "ecell4::ReactionRule::policy_type":
        Cpp_STRICT "ecell4::ReactionRule::STRICT"
        Cpp_IMPLICIT "ecell4::ReactionRule::IMPLICIT"
        Cpp_DESTROY "ecell4::ReactionRule::DESTROY"

cdef extern from "ecell4/core/ReactionRule.hpp" namespace "ecell4":
    cdef cppclass Cpp_ReactionRule "ecell4::ReactionRule":
        Cpp_ReactionRule() except +
        Cpp_ReactionRule(vector[Cpp_Species]&, vector[Cpp_Species]&)
        Cpp_ReactionRule(vector[Cpp_Species]&, vector[Cpp_Species]&, Real)
        Cpp_ReactionRule(Cpp_ReactionRule&) except +
        Real k()
        vector[Cpp_Species]& reactants()
        vector[Cpp_Species]& products()
        # multiset[Cpp_Species]& reactants()
        # multiset[Cpp_Species]& products()
        void set_k(Real)
        void add_reactant(Cpp_Species)
        void add_product(Cpp_Species)
        string as_string()
        Cpp_ReactionRulePolicyType policy()
        void set_policy(Cpp_ReactionRulePolicyType)
        Integer count(vector[Cpp_Species]) except +
        vector[Cpp_ReactionRule] generate(vector[Cpp_Species]) except +

## ReactionRule
#  a python wrapper for Cpp_ReactionRule
cdef class ReactionRule:
    cdef Cpp_ReactionRule* thisptr

cdef ReactionRule ReactionRule_from_Cpp_ReactionRule(Cpp_ReactionRule *rr)

## Cpp_Space
#  ecell4::Space
cdef extern from "ecell4/core/Space.hpp" namespace "ecell4":
    cdef cppclass Cpp_Space "ecell4::Space":
        pass

## Space
#  a python wrapper for Cpp_Space
cdef class Space:
    cdef shared_ptr[Cpp_Space]* thisptr

## Cpp_CompartmentSpaceVectorImpl
#  ecell4::CompartmentSpaceVectorImpl
cdef extern from "ecell4/core/CompartmentSpace.hpp" namespace "ecell4":
    cdef cppclass Cpp_CompartmentSpaceVectorImpl "ecell4::CompartmentSpaceVectorImpl":
        Cpp_CompartmentSpaceVectorImpl(Cpp_Real3&) except+
        Real volume()
        Integer num_molecules(Cpp_Species &sp)
        vector[Cpp_Species] list_species()
        void reset(Cpp_Real3&)
        Cpp_Real3 edge_lengths()
        void set_volume(Real)
        void add_molecules(Cpp_Species &sp, Integer num)
        void remove_molecules(Cpp_Species &sp, Integer num)

## CompartmentSpaceVectorImpl
#  a python wrapper for Cpp_CompartmentSpaceVectorImpl
cdef class CompartmentSpaceVectorImpl:
    cdef Cpp_CompartmentSpaceVectorImpl* thisptr

## Cpp_ParticleSpaceVectorImpl
#  ecell4::ParticleSpaceVectorImpl
cdef extern from "ecell4/core/ParticleSpace.hpp" namespace "ecell4":
    cdef cppclass Cpp_ParticleSpaceVectorImpl "ecell4::ParticleSpaceVectorImpl":
        Cpp_ParticleSpaceVectorImpl(Cpp_Real3&) except+
        Cpp_Real3 edge_lengths()
        Integer num_particles()
        Integer num_particles(Cpp_Species&)
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles()
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles(Cpp_Species &sp)
        bool has_particle(Cpp_ParticleID &pid)

        bool update_particle(Cpp_ParticleID, Cpp_Particle)
        pair[Cpp_ParticleID, Cpp_Particle] get_particle(Cpp_ParticleID &pid)
        void remove_particle(Cpp_ParticleID &pid)
        vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real] ] list_particles_within_radius(
                Cpp_Real3 &pos, Real &radius)
        vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real] ] list_particles_within_radius(
                Cpp_Real3 &pos, Real &radius, Cpp_ParticleID &ignore)
        vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real] ] list_particles_within_radius(
                Cpp_Real3 &pos, Real &radius, Cpp_ParticleID &ignore1, Cpp_ParticleID &ignore2)

## ParticleSpaceVectorImpl
#  a python wrapper for ParticleSpaceVectorImpl
cdef class ParticleSpaceVectorImpl:
    cdef Cpp_ParticleSpaceVectorImpl* thisptr

## Cpp_Model
#  ecell4::Model
cdef extern from "ecell4/core/Model.hpp" namespace "ecell4":
    cdef cppclass Cpp_Model "ecell4::Model":
        Cpp_Model() except +
        void add_species_attribute(Cpp_Species sp) except +
        bool has_species_attribute(Cpp_Species sp)
        void remove_species_attribute(Cpp_Species sp) except +
        void add_reaction_rule(Cpp_ReactionRule) except +
        void remove_reaction_rule(Cpp_ReactionRule) except +
        bool has_reaction_rule(Cpp_ReactionRule)
        Cpp_Species apply_species_attributes(Cpp_Species& sp)
        vector[Cpp_ReactionRule] query_reaction_rules(Cpp_Species sp)
        vector[Cpp_ReactionRule] query_reaction_rules(
            Cpp_Species sp, Cpp_Species sp)
        vector[Cpp_Species] list_species()
        # Cpp_Species create_species(string name)
        Integer num_reaction_rules()
        vector[Cpp_Species] species_attributes()
        vector[Cpp_ReactionRule] reaction_rules()

        void add_species_attributes(vector[Cpp_Species]) except +
        void add_reaction_rules(vector[Cpp_ReactionRule]) except +

        shared_ptr[Cpp_Model] expand(vector[Cpp_Species]) except +
        shared_ptr[Cpp_Model] expand(vector[Cpp_Species], Integer) except +
        shared_ptr[Cpp_Model] expand(vector[Cpp_Species], Integer, map[Cpp_Species, Integer]) except +

## Model
#  a python wrapper for Cpp_Model, but wrapped by shared_ptr
cdef class Model:
    # cdef Cpp_Model* thisptr
    # cdef shared_ptr[Cpp_Model]* thisptr
    cdef shared_ptr[Cpp_Model] thisptr

cdef Model Model_from_Cpp_Model(shared_ptr[Cpp_Model] m)

## Cpp_NetworkModel
#  ecell4::NetworkModel
cdef extern from "ecell4/core/NetworkModel.hpp" namespace "ecell4":
    cdef cppclass Cpp_NetworkModel "ecell4::NetworkModel":
        Cpp_NetworkModel() except +
        void add_species_attribute(Cpp_Species sp) except +
        bool has_species_attribute(Cpp_Species sp)
        void remove_species_attribute(Cpp_Species sp) except +
        void add_reaction_rule(Cpp_ReactionRule) except +
        void remove_reaction_rule(Cpp_ReactionRule) except +
        bool has_reaction_rule(Cpp_ReactionRule)
        Integer num_reaction_rules()
        Cpp_Species apply_species_attributes(Cpp_Species& sp)
        # Cpp_Species create_species(string name)
        vector[Cpp_Species] list_species()
        vector[Cpp_ReactionRule] query_reaction_rules(Cpp_Species sp)
        vector[Cpp_ReactionRule] query_reaction_rules(
            Cpp_Species sp, Cpp_Species sp)
        vector[Cpp_ReactionRule] reaction_rules()
        vector[Cpp_Species] species_attributes()
        void add_species_attributes(vector[Cpp_Species]) except +
        void add_reaction_rules(vector[Cpp_ReactionRule]) except +

        shared_ptr[Cpp_Model] expand(vector[Cpp_Species]) except +
        shared_ptr[Cpp_Model] expand(vector[Cpp_Species], Integer) except +
        shared_ptr[Cpp_Model] expand(vector[Cpp_Species], Integer, map[Cpp_Species, Integer]) except +

## NetworkModel
#  a python wrapper for Cpp_NetowrkModel, but wrapped by shared_ptr
cdef class NetworkModel:
    # cdef Cpp_NetworkModel* thisptr
    # cdef shared_ptr[Cpp_NetworkModel]* thisptr
    cdef shared_ptr[Cpp_NetworkModel] thisptr

cdef NetworkModel NetworkModel_from_Cpp_NetworkModel(
    shared_ptr[Cpp_NetworkModel] m)

## Cpp_NetfreeModel
#  ecell4::NetfreeModel
cdef extern from "ecell4/core/NetfreeModel.hpp" namespace "ecell4":
    cdef cppclass Cpp_NetfreeModel "ecell4::NetfreeModel":
        Cpp_NetfreeModel() except +
        void add_species_attribute(Cpp_Species sp) except +
        bool has_species_attribute(Cpp_Species sp)
        void remove_species_attribute(Cpp_Species sp) except +
        void add_reaction_rule(Cpp_ReactionRule) except +
        void remove_reaction_rule(Cpp_ReactionRule) except +
        bool has_reaction_rule(Cpp_ReactionRule)
        Integer num_reaction_rules()
        Cpp_Species apply_species_attributes(Cpp_Species& sp)
        # Cpp_Species create_species(string name)
        vector[Cpp_Species] list_species()
        vector[Cpp_ReactionRule] query_reaction_rules(Cpp_Species sp)
        vector[Cpp_ReactionRule] query_reaction_rules(
            Cpp_Species sp, Cpp_Species sp)
        vector[Cpp_ReactionRule] reaction_rules()
        vector[Cpp_Species] species_attributes()
        void add_species_attributes(vector[Cpp_Species]) except +
        void add_reaction_rules(vector[Cpp_ReactionRule]) except +

        shared_ptr[Cpp_Model] expand(vector[Cpp_Species]) except +
        shared_ptr[Cpp_Model] expand(vector[Cpp_Species], Integer) except +
        shared_ptr[Cpp_Model] expand(vector[Cpp_Species], Integer, map[Cpp_Species, Integer]) except +

        void set_effective(bool)
        bool effective()

## NetfreeModel
#  a python wrapper for Cpp_NetfreeModel, but wrapped by shared_ptr
cdef class NetfreeModel:
    # cdef Cpp_NetfreeModel* thisptr
    # cdef shared_ptr[Cpp_NetfreeModel]* thisptr
    cdef shared_ptr[Cpp_NetfreeModel] thisptr

cdef NetfreeModel NetfreeModel_from_Cpp_NetfreeModel(
    shared_ptr[Cpp_NetfreeModel] m)

# cdef shared_ptr[Cpp_Model]* Cpp_Model_from_Model(m)
cdef shared_ptr[Cpp_Model] Cpp_Model_from_Model(m)

## Cpp_Real3
#  ecell4::Real3
cdef extern from "ecell4/core/Real3.hpp" namespace "ecell4":
    cdef cppclass Cpp_Real3 "ecell4::Real3":
        Cpp_Real3() except +
        Cpp_Real3(Real, Real, Real) except +
        Cpp_Real3(Cpp_Real3 &rhs) except+
        Real& operator[](Integer)
        Cpp_Real3 operator+(Cpp_Real3, Cpp_Real3)
        Cpp_Real3 operator-(Cpp_Real3, Cpp_Real3)
        Cpp_Real3 operator/(Cpp_Real3, Real)
        Cpp_Real3 operator*(Cpp_Real3, Real)

## Real3
#  a python wrapper for Cpp_Real3
cdef class Real3:
    cdef Cpp_Real3* thisptr

cdef Real3 Real3_from_Cpp_Real3(Cpp_Real3 *p)

## Cpp_Integer3
#  ecell4::Integer3
cdef extern from "ecell4/core/Integer3.hpp" namespace "ecell4":
    cdef cppclass Cpp_Integer3 "ecell4::Integer3":
        Cpp_Integer3() except +
        Cpp_Integer3(Integer, Integer, Integer) except +
        Cpp_Integer3(Cpp_Integer3&) except +
        Integer col
        Integer row
        Integer layer
        Integer& operator[](Integer)

cdef class Integer3:
    cdef Cpp_Integer3* thisptr

cdef Integer3 Integer3_from_Cpp_Integer3(Cpp_Integer3 *g)

## Cpp_ParticleID
#  ecell4::ParticleID
cdef extern from "ecell4/core/Identifier.hpp" namespace "ecell4":
    ctypedef int lot_type
    ctypedef unsigned long long serial_type
    ctypedef pair[int, unsigned long long] value_type

    cdef cppclass Cpp_ParticleID "ecell4::ParticleID":
        Cpp_ParticleID() except+
        Cpp_ParticleID(value_type) except+
        Cpp_ParticleID(Cpp_ParticleID& rhs) except+
        Cpp_ParticleID log_add(lot_type& rhs)
        Cpp_ParticleID log_subtract(lot_type& rhs)
        Cpp_ParticleID& lot_advance(lot_type& rhs)
        Cpp_ParticleID& lot_retraace(lot_type& rhs)
        Cpp_ParticleID serial_add(serial_type& rhs)
        Cpp_ParticleID serial_subtract(serial_type& rhs)
        Cpp_ParticleID& serial_advance(serial_type& rhs)
        Cpp_ParticleID& serial_retrace(serial_type& rhs)
        # Cpp_ParticleID &operator=(Cpp_ParticleID& rhs) # XXX not yet suppoted
        bool operator==(Cpp_ParticleID& rhs)
        bool operator!=(Cpp_ParticleID& rhs)
        bool operator<(Cpp_ParticleID& rhs)
        bool operator>=(Cpp_ParticleID& rhs)
        bool operator>(Cpp_ParticleID& rhs)
        bool operator<=(Cpp_ParticleID& rhs)
        # operator value_type()
        value_type& operator() ()
        int& lot()
        unsigned long long& serial()

cdef class ParticleID:
    cdef Cpp_ParticleID* thisptr

cdef ParticleID ParticleID_from_Cpp_ParticleID(Cpp_ParticleID* p)

## Cpp_Particle
#  ecell4::Particle
cdef extern from "ecell4/core/Particle.hpp" namespace "ecell4":
    cdef cppclass Cpp_Particle "ecell4::Particle":
        Cpp_Particle() except +
        Cpp_Particle(Cpp_Species, Cpp_Real3, Real radius, Real D) except +
        Cpp_Particle(Cpp_Particle &rhs) except+
        Cpp_Real3 position()
        Real radius()
        Real D()
        Cpp_Species &species()

## Particle
#  a python wrapper for Cpp_Particle
cdef class Particle:
    cdef Cpp_Particle* thisptr

cdef Particle Particle_from_Cpp_Particle(Cpp_Particle* p)

## Cpp_Voxel
#  ecell4::Voxel
cdef extern from "ecell4/core/Voxel.hpp" namespace "ecell4":
    cdef cppclass Cpp_Voxel "ecell4::Voxel":
        Cpp_Voxel() except +
        Cpp_Voxel(Cpp_Species, Integer, Real radius, Real D) except +
        Cpp_Voxel(Cpp_Species, Integer, Real radius, Real D, string loc) except +
        Cpp_Voxel(Cpp_Voxel &rhs) except+
        Integer coordinate()
        Real D()
        Real radius()
        Cpp_Species &species()
        string loc()

## Voxel
#  a python wrapper for Cpp_Voxel
cdef class Voxel:
    cdef Cpp_Voxel* thisptr

cdef Voxel Voxel_from_Cpp_Voxel(Cpp_Voxel* p)

## Cpp_FixedIntervalNumberObserver
#  ecell4::FixedIntervalNumberObserver
cdef extern from "ecell4/core/observers.hpp" namespace "ecell4":
    cdef cppclass Cpp_Observer "ecell4::Observer":
        Real next_time()
        void reset()

    cdef cppclass Cpp_FixedIntervalNumberObserver "ecell4::FixedIntervalNumberObserver":
        Cpp_FixedIntervalNumberObserver(Real, vector[string]) except +
        Real next_time()
        Integer num_steps()
        vector[vector[Real]] data()
        vector[Cpp_Species] targets()
        void reset()
        void save(string)

    cdef cppclass Cpp_NumberObserver "ecell4::NumberObserver":
        Cpp_NumberObserver(vector[string]) except +
        Real next_time()
        Integer num_steps()
        vector[vector[Real]] data()
        vector[Cpp_Species] targets()
        void reset()
        void save(string)

    cdef cppclass Cpp_FixedIntervalHDF5Observer "ecell4::FixedIntervalHDF5Observer":
        Cpp_FixedIntervalHDF5Observer(Real, string) except +
        Real next_time()
        Integer num_steps()
        string filename()
        string filename(Integer)
        string prefix()
        void reset()

    cdef cppclass Cpp_FixedIntervalCSVObserver "ecell4::FixedIntervalCSVObserver":
        Cpp_FixedIntervalCSVObserver(Real, string) except +
        Cpp_FixedIntervalCSVObserver(Real, string, vector[string]) except +
        Real next_time()
        Integer num_steps()
        string filename()
        # void log(Cpp_Space*)
        void log(shared_ptr[Cpp_Space]&)
        void reset()
        void set_header(string&)
        void set_formatter(string&)

    cdef cppclass Cpp_CSVObserver "ecell4::CSVObserver":
        Cpp_CSVObserver(string) except +
        Cpp_CSVObserver(string, vector[string]) except +
        Real next_time()
        Integer num_steps()
        string filename()
        # void log(Cpp_Space*)
        void log(shared_ptr[Cpp_Space]&)
        void reset()
        void set_header(string&)
        void set_formatter(string&)

    cdef cppclass Cpp_FixedIntervalTrajectoryObserver "ecell4::FixedIntervalTrajectoryObserver":
        Cpp_FixedIntervalTrajectoryObserver(Real, vector[Cpp_ParticleID], bool, Real) except +
        Cpp_FixedIntervalTrajectoryObserver(Real, bool, Real) except +
        Real next_time()
        Integer num_steps()
        Integer num_tracers()
        vector[Real]& t()
        vector[vector[Cpp_Real3]] data()
        void reset()
        @staticmethod
        bool default_resolve_boundary()
        @staticmethod
        Real default_subdt()

    cdef cppclass Cpp_TimingTrajectoryObserver "ecell4::TimingTrajectoryObserver":
        Cpp_TimingTrajectoryObserver(vector[double], vector[Cpp_ParticleID], bool, Real) except +
        Cpp_TimingTrajectoryObserver(vector[double], bool, Real) except +
        Real next_time()
        Integer num_steps()
        Integer num_tracers()
        vector[Real]& t()
        vector[vector[Cpp_Real3]] data()
        void reset()
        @staticmethod
        bool default_resolve_boundary()
        @staticmethod
        Real default_subdt()

    cdef cppclass Cpp_TimingNumberObserver "ecell4::TimingNumberObserver":
        Cpp_TimingNumberObserver(vector[double], vector[string]) except +  #XXX: vector[Real]
        Real next_time()
        Integer num_steps()
        vector[vector[Real]] data()
        vector[Cpp_Species] targets()
        void reset()
        void save(string)

    cdef cppclass Cpp_TimeoutObserver "ecell4::TimeoutObserver":
        Cpp_TimeoutObserver() except +
        Cpp_TimeoutObserver(Real) except +
        Real duration()
        Real accumulation()
        Real interval()
        void reset()

    cdef cppclass Cpp_FixedIntervalTrackingObserver "ecell4::FixedIntervalTrackingObserver":
        Cpp_FixedIntervalTrackingObserver(Real, vector[Cpp_Species], bool, Real, Real) except +
        Real next_time()
        Integer num_steps()
        Integer num_tracers()
        vector[Real]& t()
        vector[vector[Cpp_Real3]] data()
        void reset()
        @staticmethod
        bool default_resolve_boundary()
        @staticmethod
        Real default_subdt()
        @staticmethod
        Real default_threshold()

## FixedIntervalNumberObserver
#  a python wrapper for Cpp_FixedIntervalNumberObserver
cdef class Observer:
    cdef shared_ptr[Cpp_Observer]* thisptr

cdef class FixedIntervalNumberObserver:
    cdef shared_ptr[Cpp_FixedIntervalNumberObserver]* thisptr

cdef class NumberObserver:
    cdef shared_ptr[Cpp_NumberObserver]* thisptr

cdef class TimingNumberObserver:
    cdef shared_ptr[Cpp_TimingNumberObserver]* thisptr

cdef class FixedIntervalHDF5Observer:
    cdef shared_ptr[Cpp_FixedIntervalHDF5Observer]* thisptr

cdef class FixedIntervalCSVObserver:
    cdef shared_ptr[Cpp_FixedIntervalCSVObserver]* thisptr

cdef class CSVObserver:
    cdef shared_ptr[Cpp_CSVObserver]* thisptr

cdef class FixedIntervalTrajectoryObserver:
    cdef shared_ptr[Cpp_FixedIntervalTrajectoryObserver]* thisptr

cdef class TimingTrajectoryObserver:
    cdef shared_ptr[Cpp_TimingTrajectoryObserver]* thisptr

cdef class TimeoutObserver:
    cdef shared_ptr[Cpp_TimeoutObserver]* thisptr

cdef class FixedIntervalTrackingObserver:
    cdef shared_ptr[Cpp_FixedIntervalTrackingObserver]* thisptr

## Cpp_Shape
#  ecell4::Shape
cdef extern from "ecell4/core/Shape.hpp" namespace "ecell4":
    cdef cppclass Cpp_Shape "ecell4::Shape":
        bool is_inside(Cpp_Real3&)
        Integer dimension()

## Cpp_Complement
#  ecell4::Complement
cdef extern from "ecell4/core/shape_operators.hpp" namespace "ecell4":
    cdef cppclass Cpp_Surface "ecell4::Surface":
        Cpp_Surface()
        Cpp_Surface(shared_ptr[Cpp_Shape]&)
        Cpp_Surface(Cpp_Surface&)
        Real is_inside(Cpp_Real3&)
        Integer dimension()

    cdef cppclass Cpp_Union "ecell4::Union":
        Cpp_Union(shared_ptr[Cpp_Shape]&, shared_ptr[Cpp_Shape]&)
        Cpp_Union(Cpp_Union&)
        Real is_inside(Cpp_Real3&)
        Integer dimension()
        Cpp_Surface surface()

    cdef cppclass Cpp_Complement "ecell4::Complement":
        Cpp_Complement(shared_ptr[Cpp_Shape]&, shared_ptr[Cpp_Shape]&)
        Cpp_Complement(Cpp_Complement&)
        Real is_inside(Cpp_Real3&)
        Integer dimension()
        Cpp_Surface surface()

    cdef cppclass Cpp_AffineTransformation "ecell4::AffineTransformation":
        Cpp_AffineTransformation()
        Cpp_AffineTransformation(shared_ptr[Cpp_Shape]&)
        Cpp_AffineTransformation(Cpp_AffineTransformation&)
        Real is_inside(Cpp_Real3&)
        Integer dimension()
        Cpp_Surface surface()
        void translate(Cpp_Real3&)
        void rescale(Cpp_Real3&)
        void xroll(Real&)
        void yroll(Real&)
        void zroll(Real&)

## Cpp_Sphere
#  ecell4::Sphere
cdef extern from "ecell4/core/Sphere.hpp" namespace "ecell4":
    cdef cppclass Cpp_Sphere "ecell4::Sphere":
        Cpp_Sphere()
        Cpp_Sphere(Cpp_Real3&, Real)
        Cpp_Sphere(Cpp_Sphere&)
        Real distance(Cpp_Real3&)
        Real is_inside(Cpp_Real3&)
        Cpp_SphericalSurface surface()
        Integer dimension()

## Cpp_SphericalSurface
#  ecell4::SphericalSurface
cdef extern from "ecell4/core/Sphere.hpp" namespace "ecell4":
    cdef cppclass Cpp_SphericalSurface "ecell4::SphericalSurface":
        Cpp_SphericalSurface()
        Cpp_SphericalSurface(Cpp_Real3&, Real)
        Cpp_SphericalSurface(Cpp_SphericalSurface&)
        Real distance(Cpp_Real3&)
        Real is_inside(Cpp_Real3&)
        Cpp_Sphere inside()
        Integer dimension()

## Cpp_Cylinder
#  ecell4::Cylinder
cdef extern from "ecell4/core/Cylinder.hpp" namespace "ecell4":
    cdef cppclass Cpp_Cylinder "ecell4::Cylinder":
        Cpp_Cylinder()
        Cpp_Cylinder(Cpp_Real3&, Real, Cpp_Real3&, Real)
        Cpp_Cylinder(Cpp_Cylinder&)
        Real distance(Cpp_Real3&)
        Real is_inside(Cpp_Real3&)
        Cpp_CylindricalSurface surface()
        Integer dimension()

## Cpp_CylindricalSurface
#  ecell4::CylindricalSurface
cdef extern from "ecell4/core/Cylinder.hpp" namespace "ecell4":
    cdef cppclass Cpp_CylindricalSurface "ecell4::CylindricalSurface":
        Cpp_CylindricalSurface()
        Cpp_CylindricalSurface(Cpp_Real3&, Real, Cpp_Real3&, Real)
        Cpp_CylindricalSurface(Cpp_CylindricalSurface&)
        Real distance(Cpp_Real3&)
        Real is_inside(Cpp_Real3&)
        Cpp_Cylinder inside()
        Integer dimension()

## Cpp_PlanarSurface
# ecell4::PlanarSurface
cdef extern from "ecell4/core/PlanarSurface.hpp" namespace "ecell4":
    cdef cppclass Cpp_PlanarSurface "ecell4::PlanarSurface":
        Cpp_PlanarSurface()
        Cpp_PlanarSurface(Cpp_Real3&, Cpp_Real3&, Cpp_Real3&)
        Cpp_PlanarSurface(Cpp_PlanarSurface)
        # Real distance(Cpp_Real3&)
        Real is_inside(Cpp_Real3&)
        Integer dimension()

## Cpp_Rod
# ecell4::Rod
cdef extern from "ecell4/core/Rod.hpp" namespace "ecell4":
    cdef cppclass Cpp_Rod "ecell4::Rod":
        Cpp_Rod()
        #Cpp_Rod(Real, Real)
        Cpp_Rod(Real, Real, Cpp_Real3&)
        Cpp_Rod(Cpp_Rod&)
        Real distance(Cpp_Real3&)
        Real is_inside(Cpp_Real3&)
        void shift(Cpp_Real3&)
        Cpp_RodSurface surface()
        Integer dimension()
        Cpp_Real3& origin()
        Real length()
        Real radius()

## Cpp_RodSurface
# ecell4::RodSurface
cdef extern from "ecell4/core/Rod.hpp" namespace "ecell4":
    cdef cppclass Cpp_RodSurface "ecell4::RodSurface":
        Cpp_RodSurface()
        #Cpp_RodSurface(Real, Real)
        Cpp_RodSurface(Real, Real, Cpp_Real3&)
        Cpp_RodSurface(Cpp_RodSurface)
        Real distance(Cpp_Real3&)
        Real is_inside(Cpp_Real3&)
        Cpp_Real3& origin()
        void shift(Cpp_Real3&)
        Cpp_Rod inside()
        Integer dimension()
        Real length()
        Real radius()

## Cpp_AABB
#  ecell4::AABB
cdef extern from "ecell4/core/AABB.hpp" namespace "ecell4":
    cdef cppclass Cpp_AABB "ecell4::AABB":
        Cpp_AABB()
        Cpp_AABB(Cpp_Real3&, Cpp_Real3&)
        Cpp_AABB(Cpp_AABB&)
        Real distance(Cpp_Real3&)
        Real is_inside(Cpp_Real3&)
        Integer dimension()
        Cpp_Real3 upper()
        Cpp_Real3 lower()
        Cpp_Surface surface()

## Cpp_MeshSurface
# ecell4::MeshSurface
cdef extern from "ecell4/core/Mesh.hpp" namespace "ecell4":
    cdef cppclass Cpp_MeshSurface "ecell4::MeshSurface":
        Cpp_MeshSurface(string, Cpp_Real3)
        Cpp_MeshSurface(Cpp_MeshSurface)
        # Real distance(Cpp_Real3&)
        Real is_inside(Cpp_Real3&)
        Integer dimension()

## Shape
#  a python wrapper for Cpp_Shape
cdef class Shape:
    cdef shared_ptr[Cpp_Shape]* thisptr

## Sphere
#  a python wrapper for Cpp_Sphere
cdef class Sphere:
    cdef shared_ptr[Cpp_Sphere]* thisptr

## SphericalSurface
#  a python wrapper for Cpp_SphericalSurface
cdef class SphericalSurface:
    cdef shared_ptr[Cpp_SphericalSurface]* thisptr

## Cylinder
#  a python wrapper for Cpp_Cylinder
cdef class Cylinder:
    cdef shared_ptr[Cpp_Cylinder]* thisptr

## CylindricalSurface
#  a python wrapper for Cpp_CylindricalSurface
cdef class CylindricalSurface:
    cdef shared_ptr[Cpp_CylindricalSurface]* thisptr

## PlanarSurface
#  a python wrapper for Cpp_PlanarSurface
cdef class PlanarSurface:
    cdef shared_ptr[Cpp_PlanarSurface]* thisptr

## Rod
# a python wrapper for Cpp_Rod
cdef class Rod:
    cdef shared_ptr[Cpp_Rod]* thisptr

## RodSurface
# a python wrapper for Cpp_RodSurface
cdef class RodSurface:
    cdef shared_ptr[Cpp_RodSurface]* thisptr


## MeshSurface
# a python wrapper for Cpp_MeshSurface
cdef class MeshSurface:
    cdef shared_ptr[Cpp_MeshSurface]* thisptr

## AABB
#  a python wrapper for Cpp_AABB
cdef class AABB:
    cdef shared_ptr[Cpp_AABB]* thisptr

## Surface
#  a python wrapper for Cpp_Surface
cdef class Surface:
    cdef shared_ptr[Cpp_Surface]* thisptr

## Union
#  a python wrapper for Cpp_Union
cdef class Union:
    cdef shared_ptr[Cpp_Union]* thisptr

## Complement
#  a python wrapper for Cpp_Complement
cdef class Complement:
    cdef shared_ptr[Cpp_Complement]* thisptr

## AffineTransformation
#  a python wrapper for Cpp_AffineTransformation
cdef class AffineTransformation:
    cdef shared_ptr[Cpp_AffineTransformation]* thisptr

cdef Sphere Sphere_from_Cpp_Sphere(Cpp_Sphere* p)
cdef SphericalSurface SphericalSurface_from_Cpp_SphericalSurface(Cpp_SphericalSurface* p)
cdef Cylinder Cylinder_from_Cpp_Cylinder(Cpp_Cylinder* p)
cdef CylindricalSurface CylindricalSurface_from_Cpp_CylindricalSurface(Cpp_CylindricalSurface* p)
cdef AABB AABB_from_Cpp_AABB(Cpp_AABB* p)

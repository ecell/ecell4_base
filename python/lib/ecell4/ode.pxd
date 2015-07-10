from libcpp.string cimport string
from libcpp cimport bool

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *


## Cpp_ODEWorld
#  ecell4::ode::ODEWorld
cdef extern from "ecell4/ode/ODEWorld.hpp" namespace "ecell4::ode":
    cdef cppclass Cpp_ODEWorld "ecell4::ode::ODEWorld":
        Cpp_ODEWorld() except +
        Cpp_ODEWorld(Cpp_Real3&) except +
        Cpp_ODEWorld(string&) except +
        # SpaceTraits
        Real& t()
        void set_t(Real&)
        void reset(Cpp_Real3&)
        Cpp_Real3 edge_lengths()
        # CompartmentSpaceTraits
        Real &volume()
        Integer num_molecules(Cpp_Species &)
        Integer num_molecules_exact(Cpp_Species &)
        vector[Cpp_Species] list_species()

        # CompartmentSpace member functions
        void set_volume(Real &)
        void add_molecules(Cpp_Species &sp, Integer &num)
        void add_molecules(Cpp_Species &sp, Integer &num, shared_ptr[Cpp_Shape])
        void remove_molecules(Cpp_Species &sp, Integer &num)
        # Optional members
        Real get_value(Cpp_Species &)
        void set_value(Cpp_Species &sp, Real &num)
        void save(string) except +
        void load(string)
        bool has_species(Cpp_Species &)
        void reserve_species(Cpp_Species &)
        void release_species(Cpp_Species &)
        void bind_to(shared_ptr[Cpp_Model])

## ODEWorld
#  a python wrapper for Cpp_ODEWorld
cdef class ODEWorld:
    cdef shared_ptr[Cpp_ODEWorld]* thisptr

cdef ODEWorld ODEWorld_from_Cpp_ODEWorld(shared_ptr[Cpp_ODEWorld] m)

## Cpp_ODESimulator
#  ecell4::ode::ODESimulator
cdef extern from "ecell4/ode/ODESimulator.hpp" namespace "ecell4::ode":
    cdef cppclass Cpp_ODESimulator "ecell4::ode::ODESimulator":
        Cpp_ODESimulator(
            shared_ptr[Cpp_Model], shared_ptr[Cpp_ODEWorld]) except +
        Cpp_ODESimulator(
            shared_ptr[Cpp_ODEWorld]) except +
        void initialize()
        Real t()
        Integer num_steps()
        Real dt()
        Real next_time()
        void step()
        bool step(Real&)
        # Optional members
        void set_t(Real&)
        void set_dt(Real &)
        vector[Cpp_ReactionRule] last_reactions()
        shared_ptr[Cpp_Model] model()
        shared_ptr[Cpp_ODEWorld] world()
        void run(Real)
        void run(Real, shared_ptr[Cpp_Observer])
        void run(Real, vector[shared_ptr[Cpp_Observer]])

## ODESimulator
#  a python wrapper for Cpp_ODESimulator
cdef class ODESimulator:
    cdef Cpp_ODESimulator *thisptr

cdef ODESimulator ODESimulator_from_Cpp_ODESimulator(Cpp_ODESimulator* s)

## Cpp_ODEFactory
#  ecell4::ode::ODEFactory
cdef extern from "ecell4/ode/ODEFactory.hpp" namespace "ecell4::ode":
    cdef cppclass Cpp_ODEFactory "ecell4::ode::ODEFactory":
        Cpp_ODEFactory() except +
        Cpp_ODEWorld* create_world(string)
        Cpp_ODEWorld* create_world(Cpp_Real3&)
        Cpp_ODEWorld* create_world(shared_ptr[Cpp_Model])
        Cpp_ODESimulator* create_simulator(shared_ptr[Cpp_Model], shared_ptr[Cpp_ODEWorld])
        Cpp_ODESimulator* create_simulator(shared_ptr[Cpp_ODEWorld])

## ODEFactory
#  a python wrapper for Cpp_ODEFactory
cdef class ODEFactory:
    cdef Cpp_ODEFactory* thisptr


## Following definitions are ODESimulator2 related.

## Cpp_ODERatelaw
cdef extern from "ecell4/ode/ODERatelaw.hpp" namespace "ecell4::ode":
    cdef cppclass Cpp_ODERatelaw "ecell4::ode::ODERatelaw":
        Cpp_ODERatelaw() except +
        bool is_available()

## ODERatelaw
cdef class ODERatelaw:
    #cdef Cpp_ODERatelaw *thisptr
    cdef shared_ptr[Cpp_ODERatelaw] *thisptr

## Cpp_ODERatelawMassAction
cdef extern from "ecell4/ode/ODERatelaw.hpp" namespace "ecell4::ode":
    cdef cppclass Cpp_ODERatelawMassAction "ecell4::ode::ODERatelawMassAction":
        Cpp_ODERatelawMassAction(Real) except +
        bool is_available()
        void set_k(Real)
        Real get_k()

cdef class ODERatelawMassAction:
    #cdef Cpp_ODERatelawMassAction *thisptr
    cdef shared_ptr[Cpp_ODERatelawMassAction] *thisptr

ctypedef void* Python_CallbackFunctype
ctypedef double (*Stepladder_Functype)(
    Python_CallbackFunctype pyfunc, vector[Real], vector[Real], 
    Real volume, Real t, Cpp_ODEReactionRule *)

cdef extern from "ecell4/ode/ODERatelaw.hpp" namespace "ecell4::ode":
    cdef cppclass Cpp_ODERatelawCythonCallback " ecell4::ode::ODERatelawCythonCallback":
        Cpp_ODERatelawCythonCallback() except+
        Cpp_ODERatelawCythonCallback(Stepladder_Functype, Python_CallbackFunctype) except+
        bool is_available()
        void set_callback_pyfunc(Python_CallbackFunctype)

cdef class ODERatelawCallback:
    cdef shared_ptr[Cpp_ODERatelawCythonCallback] *thisptr
    cdef object pyfunc

## Cpp_ODEReactionRule
cdef extern from "ecell4/ode/ODEReactionRule.hpp" namespace "ecell4::ode":
    cdef cppclass Cpp_ODEReactionRule "ecell4::ode::ODEReactionRule":
        Cpp_ODEReactionRule() except +
        Cpp_ODEReactionRule(Cpp_ReactionRule) except +
        Cpp_ODEReactionRule(Cpp_ODEReactionRule) except +
        Real k()
        void set_k(Real)
        vector[Cpp_Species] reactants()
        vector[Cpp_Species] products()
        vector[Real] reactants_coefficients()
        vector[Real] products_coefficients()

        void add_reactant(Cpp_Species, Real)
        void add_product(Cpp_Species, Real)
        void set_reactant_coefficient(int, Real)
        void set_product_coefficient(int, Real)

        void set_ratelaw(shared_ptr[Cpp_ODERatelaw])
        void set_ratelaw(shared_ptr[Cpp_ODERatelawMassAction])
        shared_ptr[Cpp_ODERatelaw] get_ratelaw()
        bool has_ratelaw()
        bool is_massaction()

cdef class ODEReactionRule:
    cdef Cpp_ODEReactionRule *thisptr
        
## Cpp_ODENetworkModel
cdef extern from "ecell4/ode/ODENetworkModel.hpp" namespace "ecell4::ode":
    cdef cppclass Cpp_ODENetworkModel "ecell4::ode::ODENetworkModel":
        Cpp_ODENetworkModel() except +
        Cpp_ODENetworkModel( shared_ptr[Cpp_NetworkModel] ) except +
        void update_model()
        bool has_model()
        vector[Cpp_ODEReactionRule] ode_reaction_rules()
        Integer num_reaction_rules()
        void dump_reactions()
        void add_reaction_rule(Cpp_ODEReactionRule)

cdef class ODENetworkModel:
    cdef shared_ptr[Cpp_ODENetworkModel] *thisptr

## Cpp_ODESimulator2
cdef extern from "ecell4/ode/ODESimulator2.hpp" namespace "ecell4::ode":
    cdef cppclass Cpp_ODESimulator2 "ecell4::ode::ODESimulator2":
        Cpp_ODESimulator2(shared_ptr[Cpp_ODENetworkModel], shared_ptr[Cpp_ODEWorld]) except+
        void initialize()
        void step()
        bool step(Real)
        Real next_time()
        Real t()
        void set_t(Real)
        Real dt()
        void set_dt(Real)
        Integer num_steps()

cdef class ODESimulator2:
    cdef Cpp_ODESimulator2 *thisptr

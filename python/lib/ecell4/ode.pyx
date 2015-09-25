import collections
from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address
from libcpp.string cimport string
from libcpp.vector cimport vector

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *

from cpython cimport PyObject, Py_XINCREF, Py_XDECREF

## ODEWorld
#  a python wrapper for Cpp_ODEWorld
cdef class ODEWorld:

    def __cinit__(self, edge_lengths = None):
        cdef string filename

        if edge_lengths is None:
            self.thisptr = new shared_ptr[Cpp_ODEWorld](new Cpp_ODEWorld())
        elif isinstance(edge_lengths, Real3):
            self.thisptr = new shared_ptr[Cpp_ODEWorld](
                new Cpp_ODEWorld(deref((<Real3>edge_lengths).thisptr)))
        else:
            filename = tostring(edge_lengths)
            self.thisptr = new shared_ptr[Cpp_ODEWorld](new Cpp_ODEWorld(filename))

    def __dealloc__(self):
        # XXX: Here, we release shared pointer,
        #      and if reference count to the ODEWorld object become zero,
        #      it will be released automatically.
        del self.thisptr

    def set_t(self, Real t):
        self.thisptr.get().set_t(t)

    def t(self):
        return self.thisptr.get().t()

    def edge_lengths(self):
        cdef Cpp_Real3 lengths = self.thisptr.get().edge_lengths()
        return Real3_from_Cpp_Real3(address(lengths))

    def volume(self):
        return self.thisptr.get().volume()

    def num_molecules(self, Species sp):
        return self.thisptr.get().num_molecules(deref(sp.thisptr))

    def num_molecules_exact(self, Species sp):
        return self.thisptr.get().num_molecules_exact(deref(sp.thisptr))

    def list_species(self):
        cdef vector[Cpp_Species] raw_list_species = self.thisptr.get().list_species()
        retval = []
        cdef vector[Cpp_Species].iterator it = raw_list_species.begin()
        while it != raw_list_species.end():
            retval.append(
                Species_from_Cpp_Species(<Cpp_Species*> (address(deref(it)))))
            inc(it)
        return retval

    def set_volume(self, Real vol):
        self.thisptr.get().set_volume(vol)

    def add_molecules(self, Species sp, Integer num, shape=None):
        if shape is None:
            self.thisptr.get().add_molecules(deref(sp.thisptr), num)
        else:
            self.thisptr.get().add_molecules(
                deref(sp.thisptr), num, deref((<Shape>(shape.as_base())).thisptr))

    def remove_molecules(self, Species sp, Integer num):
        self.thisptr.get().remove_molecules(deref(sp.thisptr), num)

    def get_value(self, Species sp):
        return self.thisptr.get().get_value(deref(sp.thisptr))

    def set_value(self, Species sp, Real num):
        self.thisptr.get().set_value(deref(sp.thisptr), num)

    def save(self, filename):
        self.thisptr.get().save(tostring(filename))

    def load(self, string filename):
        self.thisptr.get().load(tostring(filename))

    def has_species(self, Species sp):
        return self.thisptr.get().has_species(deref(sp.thisptr))

    def reserve_species(self, Species sp):
        self.thisptr.get().reserve_species(deref(sp.thisptr))

    def release_species(self, Species sp):
        self.thisptr.get().release_species(deref(sp.thisptr))

    def bind_to(self, m):
        if isinstance(m, ODENetworkModel):
            self.thisptr.get().bind_to(deref((<ODENetworkModel>m).thisptr))
        else:
            self.thisptr.get().bind_to(deref(Cpp_Model_from_Model(m)))

    def as_base(self):
        retval = Space()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Space](
            <shared_ptr[Cpp_Space]>deref(self.thisptr))
        return retval

cdef ODEWorld ODEWorld_from_Cpp_ODEWorld(
    shared_ptr[Cpp_ODEWorld] w):
    r = ODEWorld(Real3(1, 1, 1))
    r.thisptr.swap(w)
    return r

cdef class ODERatelaw:
    # Abstract ODERatelaw Type.
    def __cinit__(self):
        self.thisptr = new shared_ptr[Cpp_ODERatelaw](
                <Cpp_ODERatelaw*>(new Cpp_ODERatelawMassAction(0.0)) )  # Dummy

    def __dealloc__(self):
        del self.thisptr

    def as_base(self):
        return self

cdef class ODERatelawMassAction:
    def __cinit__(self, Real k):
        self.thisptr = new shared_ptr[Cpp_ODERatelawMassAction](
                <Cpp_ODERatelawMassAction*>(new Cpp_ODERatelawMassAction(k)))

    def __dealloc__(self):
        del self.thisptr

    def is_available(self):
        return self.get().is_available
    def set_k(self, Real k):
        self.get().thisptr.set_k(k)
    def get_k(self):
        return self.get().thisptr.get_k()
    def as_base(self):
        base_type = ODERatelaw()
        del base_type.thisptr
        base_type.thisptr = new shared_ptr[Cpp_ODERatelaw](
                <shared_ptr[Cpp_ODERatelaw]>(deref(self.thisptr))
                )
        return base_type


cdef double indirect_function(
    void *func, vector[Real] reactants, vector[Real] products, 
    Real volume, Real t, Cpp_ODEReactionRule *rr):
    py_reactants = []
    cdef vector[Real].iterator it1 = reactants.begin()
    while it1 != reactants.end():
        py_reactants.append(deref(it1))
        inc(it1)
    py_products = []
    cdef vector[Real].iterator it2 = products.begin()
    while it2 != products.end():
        py_products.append(deref(it2))
        inc(it2)
    return (<object>func)(
            py_reactants, py_products, volume, t, ODEReactionRule_from_Cpp_ODEReactionRule(rr))

cdef void inc_ref(void* func):
    Py_XINCREF(<PyObject*>func)

cdef void dec_ref(void* func):
    Py_XDECREF(<PyObject*>func)

cdef class ODERatelawCallback:
    def __cinit__(self, pyfunc):
        self.thisptr = new shared_ptr[Cpp_ODERatelawCythonCallback](
            <Cpp_ODERatelawCythonCallback*>(new Cpp_ODERatelawCythonCallback(
                <Stepladder_Functype>indirect_function, <void*>pyfunc, 
                <OperateRef_Functype>inc_ref, <OperateRef_Functype>dec_ref)))
        self.pyfunc = pyfunc

    def __dealloc__(self):
        del self.thisptr

    def set_callback(self, pyfunc):
        self.thisptr.get().set_callback_pyfunc(<Python_CallbackFunctype>pyfunc)

    def as_base(self):
        retval = ODERatelaw()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_ODERatelaw](
            <shared_ptr[Cpp_ODERatelaw]>deref(self.thisptr))
        return retval

cdef class ODEReactionRule:
    def __cinit__(self):
        self.thisptr = new Cpp_ODEReactionRule()
        self.ratelaw = None

    def __dealloc__(self):
        del self.thisptr

    def k(self):
        return self.thisptr.k()
    def set_k(self, Real k):
        self.thisptr.set_k(k)

    def add_reactant(self, Species sp, coeff=None):
        if coeff is not None:
            self.thisptr.add_reactant(deref(sp.thisptr), coeff)
        else:
            self.thisptr.add_reactant(deref(sp.thisptr))

    def add_product(self, Species sp, coeff=None):
        if coeff is not None:
            self.thisptr.add_product(deref(sp.thisptr), coeff)
        else:
            self.thisptr.add_product(deref(sp.thisptr))

    def set_reactant_coefficient(self, Integer index, Real coeff):
        self.thisptr.set_reactant_coefficient(index, coeff)
    def set_product_coefficient(self, Integer index, Real coeff):
        self.thisptr.set_product_coefficient(index, coeff)

    def set_ratelaw(self, ratelaw_obj):
        self.ratelaw = ratelaw_obj
        self.thisptr.set_ratelaw(deref( (<ODERatelaw>(ratelaw_obj.as_base())).thisptr )) 

    def set_ratelaw_massaction(self, ODERatelawMassAction ratelaw_obj):
        self.ratelaw = ratelaw_obj
        self.thisptr.set_ratelaw( deref(ratelaw_obj.thisptr) )

    def has_ratelaw(self):
        return self.thisptr.has_ratelaw()

    def is_massaction(self):
        return self.thisptr.is_massaction()

    def reactants(self):
        cdef vector[Cpp_Species] cpp_reactants = self.thisptr.reactants()
        retval = []
        cdef vector[Cpp_Species].iterator it = cpp_reactants.begin()
        while it != cpp_reactants.end():
            retval.append(
                    Species_from_Cpp_Species(<Cpp_Species*>address(deref(it))))
            inc(it)
        return retval

    def reactants_coefficients(self):
        cdef vector[Real] coefficients = self.thisptr.reactants_coefficients()
        retval = []
        cdef vector[Real].iterator it = coefficients.begin()
        while it != coefficients.end():
            retval.append( deref(it) )
            inc(it)
        return retval

    def products(self):
        cdef vector[Cpp_Species] cpp_products = self.thisptr.products()
        retval = []
        cdef vector[Cpp_Species].iterator it = cpp_products.begin()
        while it != cpp_products.end():
            retval.append(
                    Species_from_Cpp_Species(<Cpp_Species*>address(deref(it))))
            inc(it)
        return retval

    def products_coefficients(self):
        cdef vector[Real] coefficients = self.thisptr.products_coefficients()
        retval = []
        cdef vector[Real].iterator it = coefficients.begin()
        while it != coefficients.end():
            retval.append( deref(it) )
            inc(it)
        return retval

    def as_string(self):
        reactants = self.reactants()
        reactants_coeff = self.reactants_coefficients()
        products = self.products()
        products_coeff = self.products_coefficients()
        leftside = ""
        rightside = ""
        retval = ""
        first = True
        for (sp, coeff) in zip(reactants, reactants_coeff):
            s = "{0}({1})".format(coeff, sp.serial())
            if first == True:
                leftside = s
                first = False
            else:
                leftside = "{} + {}".format(leftside, s)
        first = True
        for (sp, coeff) in zip(products, products_coeff):
            s = "{0}({1})".format(coeff, sp.serial())
            if first == True:
                rightside += s
                first = False
            else:
                rightside = "{} + {}".format(retval, s)
        s = ""
        if self.has_ratelaw():
            s = "HAVE"
        else:
            s = "DON'T HAVE"
        if self.is_massaction():
            k_desc = "k = {:f}\t {} Ratelaw".format(self.k(), s)
        else:
            k_desc = "\t {} Ratelaw".format(s)
        retval = "{} ---> {}\t{}".format(leftside, rightside, k_desc)
        return retval

cdef ODEReactionRule ODEReactionRule_from_Cpp_ODEReactionRule(Cpp_ODEReactionRule *s):
    cdef Cpp_ODEReactionRule *new_obj = new Cpp_ODEReactionRule(deref(s))
    ret = ODEReactionRule()
    del ret.thisptr
    ret.thisptr = new_obj
    return ret

cdef class ODENetworkModel:
    #def __cinit__(self):
    #    self.thisptr = new shared_ptr[Cpp_ODENetworkModel](
    #        <Cpp_ODENetworkModel*>(new Cpp_ODENetworkModel()) )
    def __cinit__(self, NetworkModel m = None):
        if m == None:
            self.thisptr = new shared_ptr[Cpp_ODENetworkModel](
                <Cpp_ODENetworkModel*>(new Cpp_ODENetworkModel()) )
        else:
            self.thisptr = new shared_ptr[Cpp_ODENetworkModel](
                (<Cpp_ODENetworkModel*>(new Cpp_ODENetworkModel( deref(m.thisptr) ) )) )

    def __dealloc__(self):
        del self.thisptr
    def update_model(self):
        self.thisptr.get().update_model()
    def has_network_model(self):
        return self.thisptr.get().has_network_model()
    def ode_reaction_rules(self):
        cdef vector[Cpp_ODEReactionRule] cpp_rules = self.thisptr.get().ode_reaction_rules()
        retval = []
        cdef vector[Cpp_ODEReactionRule].iterator it = cpp_rules.begin()
        while it != cpp_rules.end():
            retval.append( ODEReactionRule_from_Cpp_ODEReactionRule(address(deref(it))) )
            inc(it)
        return retval
    def num_reaction_rules(self):
        return self.thisptr.get().num_reaction_rules()

    def add_reaction_rule(self, rr):
        if isinstance(rr, ODEReactionRule):
            self.thisptr.get().add_reaction_rule(deref((<ODEReactionRule>rr).thisptr))
        elif isinstance(rr, ReactionRule):
            self.thisptr.get().add_reaction_rule(deref((<ReactionRule>rr).thisptr))
        else:
            raise ValueError("invalid argument {}".format(repr(rr)))

    def list_species(self):
        cdef vector[Cpp_Species] species = self.thisptr.get().list_species()
        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(Species_from_Cpp_Species(
                <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

cdef ODENetworkModel ODENetworkModel_from_Cpp_ODENetworkModel(
    shared_ptr[Cpp_ODENetworkModel] m):
    r = ODENetworkModel()
    r.thisptr.swap(m)
    return r

# ODESolverType:
(
    RUNGE_KUTTA_CASH_KARP54,
    ROSENBROCK4,
    EULER,
) = (0, 1, 2)

cdef Cpp_ODESolverType translate_solver_type(solvertype_constant):
    if solvertype_constant == RUNGE_KUTTA_CASH_KARP54:
        return RUNGE_KUTTA_CASH_KARP54
    elif solvertype_constant == ROSENBROCK4:
        return Cpp_ROSENBROCK4
    elif solvertype_constant == EULER:
        return Cpp_EULER
    else:
        raise ValueError(
            "invalid solver type was given [{0}]".format(repr(solvertype_constant)))

cdef class ODESimulator:

    def __cinit__(self, arg1, arg2 = None, arg3 = None):
        if arg2 is None:
            self.thisptr = new Cpp_ODESimulator(deref((<ODEWorld>arg1).thisptr))
        elif arg3 is None:
            if isinstance(arg2, ODEWorld):
                if isinstance(arg1, ODENetworkModel):
                    self.thisptr = new Cpp_ODESimulator(
                        deref((<ODENetworkModel>arg1).thisptr),
                        deref((<ODEWorld>arg2).thisptr))
                elif isinstance(arg1, NetworkModel):
                    self.thisptr = new Cpp_ODESimulator(
                        deref((<NetworkModel>arg1).thisptr),
                        deref((<ODEWorld>arg2).thisptr))
                else:
                    raise ValueError(
                        "An invalid value [{}] for the first argument.".format(repr(arg1))
                        + " NetworkModel or ODENetworkModel is needed.")
            elif isinstance(arg1, ODEWorld):
                # arg2: ODESolverType
                self.thisptr = new Cpp_ODESimulator(
                    deref((<ODEWorld>arg1).thisptr), translate_solver_type(arg2))
            else:
                raise ValueError(
                    "An invalid value [{}] for the first argument.".format(repr(arg1))
                    + " ODEWorld is needed.")
        else:
            if not isinstance(arg2, ODEWorld):
                raise ValueError(
                    "An invalid argument [{}] for the second argument.".format(repr(arg2))
                    + " ODEWorld is needed.")

            cpp_solvertype = translate_solver_type(arg3)

            if isinstance(arg1, ODENetworkModel):
                self.thisptr = new Cpp_ODESimulator(
                    deref((<ODENetworkModel>arg1).thisptr),
                    deref((<ODEWorld>arg2).thisptr),
                    cpp_solvertype)
            elif isinstance(arg1, NetworkModel):
                self.thisptr = new Cpp_ODESimulator(
                    deref((<NetworkModel>arg1).thisptr),
                    deref((<ODEWorld>arg2).thisptr),
                    cpp_solvertype)
            else:
                raise ValueError(
                    "An invalid value [{}] for the first argument.".format(repr(arg1))
                    + " NetworkModel or ODENetworkModel is needed.")

    # def __cinit__(self, m, ODEWorld w, solvertype = None):
    #     if solvertype is None:
    #         if isinstance(m, ODENetworkModel):
    #             self.thisptr = new Cpp_ODESimulator(
    #                 deref((<ODENetworkModel>m).thisptr), deref(w.thisptr))
    #         elif isinstance(m, NetworkModel):
    #             self.thisptr = new Cpp_ODESimulator(
    #                 deref((<NetworkModel>m).thisptr), deref(w.thisptr))
    #         else:
    #             raise ValueError(
    #                 "invalid argument {}.".format(repr(m))
    #                 + " NetworkModel or ODENetworkModel is needed.")
    #     else:
    #         cpp_solvertype = translate_solver_type(solvertype)

    #         if isinstance(m, ODENetworkModel):
    #             self.thisptr = new Cpp_ODESimulator(
    #                 deref((<ODENetworkModel>m).thisptr), deref(w.thisptr), cpp_solvertype)
    #         elif isinstance(m, NetworkModel):
    #             self.thisptr = new Cpp_ODESimulator(
    #                 deref((<NetworkModel>m).thisptr), deref(w.thisptr), cpp_solvertype)
    #         else:
    #             raise ValueError(
    #                 "invalid argument {}.".format(repr(m))
    #                 + " NetworkModel or ODENetworkModel is needed.")

    def __dealloc__(self):
        del self.thisptr
    def initialize(self):
        self.thisptr.initialize()
    def step(self, upto = None):
        if upto == None:
            self.thirptr.step()
        else:
            return self.thisptr.step(upto)

    def next_time(self):
        return self.thisptr.next_time()
    def t(self):
        return self.thisptr.t()
    def set_t(self, Real t_new):
        self.thisptr.set_t(t_new)
    def dt(self):
        return self.thisptr.dt()
    def set_dt(self, dt_new):
        self.thisptr.set_dt(dt_new)
    def num_steps(self):
        return self.thisptr.num_steps()
    def absolute_tolerance(self):
        return self.thisptr.absolute_tolerance()
    def set_absolute_tolerance(self, Real abs_tol):
        self.thisptr.set_absolute_tolerance(abs_tol)
    def relative_tolerance(self):
        return self.thisptr.relative_tolerance()
    def set_relative_tolerance(self, Real rel_tol):
        self.thisptr.set_relative_tolerance(rel_tol)
    def run(self, Real duration, observers=None):
        cdef vector[shared_ptr[Cpp_Observer]] tmp

        if observers is None:
            self.thisptr.run(duration)
        elif isinstance(observers, collections.Iterable):
            for obs in observers:
                tmp.push_back(deref((<Observer>(obs.as_base())).thisptr))
            self.thisptr.run(duration, tmp)
        else:
            self.thisptr.run(duration,
                deref((<Observer>(observers.as_base())).thisptr))

cdef ODESimulator ODESimulator_from_Cpp_ODESimulator(Cpp_ODESimulator* s):
    r = ODESimulator(
        ODENetworkModel_from_Cpp_ODENetworkModel(s.model()),
        ODEWorld_from_Cpp_ODEWorld(s.world()))
    del r.thisptr
    r.thisptr = s
    return r

## ODEFactory
#  a python wrapper for Cpp_ODEFactory
cdef class ODEFactory:

    def __cinit__(self, solvertype = None, dt = None):
        if solvertype is None:
            self.thisptr = new Cpp_ODEFactory()
        elif dt is None:
            self.thisptr = new Cpp_ODEFactory(translate_solver_type(solvertype))
        else:
            self.thisptr = new Cpp_ODEFactory(translate_solver_type(solvertype), dt)

    def __dealloc__(self):
        del self.thisptr

    def create_world(self, arg1=None):
        if arg1 is None:
            return ODEWorld_from_Cpp_ODEWorld(
                shared_ptr[Cpp_ODEWorld](self.thisptr.create_world()))
        elif isinstance(arg1, Real3):
            return ODEWorld_from_Cpp_ODEWorld(
                shared_ptr[Cpp_ODEWorld](
                    self.thisptr.create_world(deref((<Real3>arg1).thisptr))))
        elif isinstance(arg1, str):
            return ODEWorld_from_Cpp_ODEWorld(
                shared_ptr[Cpp_ODEWorld](self.thisptr.create_world(<string>(arg1))))
        raise ValueError("invalid argument")

    # def create_simulator(self, arg1, ODEWorld arg2=None):
    #     if arg2 is None:
    #         return ODESimulator_from_Cpp_ODESimulator(
    #             self.thisptr.create_simulator(deref((<ODEWorld>arg1).thisptr)))
    #     else:
    #         return ODESimulator_from_Cpp_ODESimulator(
    #             self.thisptr.create_simulator(
    #                 deref((<ODENetworkModel>arg1).thisptr), deref(arg2.thisptr)))

    def create_simulator(self, arg1, arg2 = None):
        if arg2 is None:
            if isinstance(arg1, ODEWorld):
                return ODESimulator_from_Cpp_ODESimulator(
                    self.thisptr.create_simulator(
                        deref((<ODEWorld>arg1).thisptr)))
            else:
                raise ValueError(
                    "invalid argument {}.".format(repr(arg1))
                    + " ODEWorld is needed.")
        else:
            if isinstance(arg1, ODENetworkModel):
                return ODESimulator_from_Cpp_ODESimulator(
                    self.thisptr.create_simulator(
                        deref((<ODENetworkModel>arg1).thisptr),
                        deref((<ODEWorld>arg2).thisptr)))
            elif isinstance(arg1, NetworkModel):
                return ODESimulator_from_Cpp_ODESimulator(
                    self.thisptr.create_simulator(
                        deref((<NetworkModel>arg1).thisptr),
                        deref((<ODEWorld>arg2).thisptr)))
            else:
                raise ValueError(
                    "invalid argument {}.".format(repr(arg1))
                    + " NetworkModel or ODENetworkModel is needed.")

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
    """A class representing the World for ODE simulations.

    ODEWorld(edge_lengths=None)

    """

    def __init__(self, edge_lengths = None):
        """Constructor.

        Args:
            edge_lengths (Real3, optional): A size of the World.

        """
        pass  # XXX: Only used for doc string

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
        """set_t(t)

        Set the current time."""
        self.thisptr.get().set_t(t)

    def t(self):
        """Return the current time."""
        return self.thisptr.get().t()

    def edge_lengths(self):
        """edge_lengths() -> Real3

        Return edge lengths for the space."""
        cdef Cpp_Real3 lengths = self.thisptr.get().edge_lengths()
        return Real3_from_Cpp_Real3(address(lengths))

    def set_volume(self, Real vol):
        """set_volume(volume)

        Set a volume."""
        self.thisptr.get().set_volume(vol)

    def volume(self):
        """Return a volume."""
        return self.thisptr.get().volume()

    def num_molecules(self, Species sp):
        """num_molecules(sp) -> Integer

        Return the number of molecules. A value is rounded to an integer.
        See set_value also.

        Args:
            sp (Species, optional): a species whose molecules you count

        Returns:
            Integer: the number of molecules (of a given species)

        """
        return self.thisptr.get().num_molecules(deref(sp.thisptr))

    def num_molecules_exact(self, Species sp):
        """num_molecules_exact(sp) -> Integer

        Return the number of molecules of a given species.
        A value is rounded to an integer. See get_value_exact also.

        Args:
            sp (Species): a species whose molecules you count

        Returns:
            Integer: the number of molecules of a given species

        """
        return self.thisptr.get().num_molecules_exact(deref(sp.thisptr))

    def list_species(self):
        """Return a list of species."""
        cdef vector[Cpp_Species] raw_list_species = self.thisptr.get().list_species()
        retval = []
        cdef vector[Cpp_Species].iterator it = raw_list_species.begin()
        while it != raw_list_species.end():
            retval.append(
                Species_from_Cpp_Species(<Cpp_Species*> (address(deref(it)))))
            inc(it)
        return retval

    def add_molecules(self, Species sp, Integer num, shape=None):
        """add_molecules(sp, num, shape=None)

        Add some molecules.

        Args:
            sp (Species): a species of molecules to add
            num (Integer): the number of molecules to add
            shape (Shape, optional): a shape to add molecules on
                [not supported yet]

        """
        if shape is None:
            self.thisptr.get().add_molecules(deref(sp.thisptr), num)
        else:
            self.thisptr.get().add_molecules(
                deref(sp.thisptr), num, deref((<Shape>(shape.as_base())).thisptr))

    def remove_molecules(self, Species sp, Integer num):
        """remove_molecules(sp, num)

        Remove molecules

        Args:
            sp (Species): a species whose molecules to remove
            num (Integer): a number of molecules to be removed

        """
        self.thisptr.get().remove_molecules(deref(sp.thisptr), num)

    def get_value(self, Species sp):
        """get_value(sp) -> Real

        Return the value matched to a given species.

        Args:
            sp (Species): a pattern whose value you get

        Returns:
            Real: the value matched to a given species

        """
        return self.thisptr.get().get_value(deref(sp.thisptr))

    def get_value_exact(self, Species sp):
        """get_value_exact(sp) -> Real

        Return the value connected to a given species.

        Args:
            sp (Species): a species whose value you get

        Returns:
            Real: the value connected to a given species

        """
        return self.thisptr.get().get_value(deref(sp.thisptr))

    def set_value(self, Species sp, Real value):
        """set_value(sp, value)

        Set the value of the given species.

        Args:
            sp (Species): a species whose value you set
            value (Real): a value set

        """
        self.thisptr.get().set_value(deref(sp.thisptr), value)

    def save(self, filename):
        """save(filename)

        Save the current state to a HDF5 file.

        Args:
            filename (str): a file name to be saved.

        """
        self.thisptr.get().save(tostring(filename))

    def load(self, string filename):
        """load(filename)

        Load a HDF5 file to the current state.

        Args:
            filename (str): a file name to be loaded.

        """
        self.thisptr.get().load(tostring(filename))

    def has_species(self, Species sp):
        """has_species(sp) -> bool

        Check if the given species is belonging to this.

        Args:
            sp (Species): a species to be checked.

        Returns:
            bool: True if the given species is contained.

        """
        return self.thisptr.get().has_species(deref(sp.thisptr))

    def reserve_species(self, Species sp):
        """reserve_species(sp)

        Reserve a value for the given species. Use set_value.

        Args:
            sp (Species): a species to be reserved.

        """
        self.thisptr.get().reserve_species(deref(sp.thisptr))

    def release_species(self, Species sp):
        """release_species(sp)

        Release a value for the given species.
        This function is mainly for developers.

        Args:
            sp (Species): a species to be released.

        """
        self.thisptr.get().release_species(deref(sp.thisptr))

    def bind_to(self, m):
        """bind_to(m)

        Bind a model.

        Args:
            m (ODENetworkModel or NetworkModel): a model to be bound

        """
        if isinstance(m, ODENetworkModel):
            self.thisptr.get().bind_to(deref((<ODENetworkModel>m).thisptr))
        else:
            self.thisptr.get().bind_to(deref(Cpp_Model_from_Model(m)))

    def as_base(self):
        """Return self as a base class. Only for developmental use."""
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
    """An abstract base class for ratelaws bound to ODEReactionRule.

    ODERatelaw()

    """

    def __init__(self):
        """Constructor."""
        pass

    def __cinit__(self):
        self.thisptr = new shared_ptr[Cpp_ODERatelaw](
                <Cpp_ODERatelaw*>(new Cpp_ODERatelawMassAction(0.0)))  # Dummy

    def __dealloc__(self):
        del self.thisptr

    def as_base(self):
        """Return self as a base class. Only for developmental use."""
        return self

cdef class ODERatelawMassAction:
    """A class for mass action ratelaws.

    ODERatelawMassAction(Real k)

    """

    def __init__(self, Real k):
        """Constructor.

        Args:
            k (Real): A kinetic rate for the mass action.

        """
        pass

    def __cinit__(self, Real k):
        self.thisptr = new shared_ptr[Cpp_ODERatelawMassAction](
                <Cpp_ODERatelawMassAction*>(new Cpp_ODERatelawMassAction(k)))

    def __dealloc__(self):
        del self.thisptr

    def is_available(self):
        """Check if this ratelaw is available or not. Return True always."""
        return self.get().is_available

    def set_k(self, Real k):
        """set_k(k)

        Set a kinetic rate constant.

        Args:
            k (float): A kinetic rate constant.

        """
        self.get().thisptr.set_k(k)

    def get_k(self):
        """Return the kinetic rate constant as a float value."""
        return self.get().thisptr.get_k()

    def as_base(self):
        """Return self as a base class. Only for developmental use."""
        base_type = ODERatelaw()
        del base_type.thisptr
        base_type.thisptr = new shared_ptr[Cpp_ODERatelaw](
                <shared_ptr[Cpp_ODERatelaw]>(deref(self.thisptr)))
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
            py_reactants, py_products, volume, t,
            ODEReactionRule_from_Cpp_ODEReactionRule(rr))

cdef void inc_ref(void* func):
    Py_XINCREF(<PyObject*>func)

cdef void dec_ref(void* func):
    Py_XDECREF(<PyObject*>func)

cdef class ODERatelawCallback:
    """A class for general ratelaws with a callback.

    ODERatelawCallback(pyfunc)

    """

    def __init__(self, pyfunc):
        """Constructor.

        Args:
            pyfunc: A Python function for the callback.
                See set_callback function of this class for details.

        """
        pass

    def __cinit__(self, pyfunc):
        self.thisptr = new shared_ptr[Cpp_ODERatelawCythonCallback](
            <Cpp_ODERatelawCythonCallback*>(new Cpp_ODERatelawCythonCallback(
                <Stepladder_Functype>indirect_function, <void*>pyfunc,
                <OperateRef_Functype>inc_ref, <OperateRef_Functype>dec_ref)))
        self.pyfunc = pyfunc

    def __dealloc__(self):
        del self.thisptr

    def set_callback(self, pyfunc):
        """set_callback(pyfunc)

        Args:
            pyfunc: A Python function for the callback
                The function must accept five arguments, and return a velocity.
                The number of reactants, the number of products, a volume,
                the current time, and a ODEReactionRule are given as the
                arguments in this order.

        Example:
            The following callback represents a simple Michaelis-Menten-like
            equation:

            >>> rl = ODERatelawCallback()
            >>> rl.set_callback(lambda r, p, v, t, rr: 2.0 * r[0] * r[1] / (1.0 + r[1]))

            Here, we expect that the first reactant is an enzyme,
            and that the second one is a substrate.

        """
        self.thisptr.get().set_callback_pyfunc(<Python_CallbackFunctype>pyfunc)

    def as_base(self):
        """Return self as a base class. Only for developmental use."""
        retval = ODERatelaw()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_ODERatelaw](
            <shared_ptr[Cpp_ODERatelaw]>deref(self.thisptr))
        return retval

cdef class ODEReactionRule:
    """A class representing a reaction rule between ``Species``, which accepts at most
    one rate law to calculate the flux.

    ODEReactionRule()

    """

    def __init__(self):
        """Constructor."""
        pass

    def __cinit__(self):
        self.thisptr = new Cpp_ODEReactionRule()
        self.ratelaw = None

    def __dealloc__(self):
        del self.thisptr

    def k(self):
        """Return the kinetic rate constant as a float value."""
        return self.thisptr.k()

    def set_k(self, Real k):
        """set_k(k)

        Set a kinetic rate constant.

        Args:
            k (float): A kinetic rate constant.

        """
        self.thisptr.set_k(k)

    def add_reactant(self, Species sp, coeff=None):
        """add_reactant(sp, coeff=None)

        Append a reactant to the end.

        Args:
            sp (Species): A new reactant.
            coeff (Integer): A stoichiometry coefficient.

        """
        if coeff is not None:
            self.thisptr.add_reactant(deref(sp.thisptr), coeff)
        else:
            self.thisptr.add_reactant(deref(sp.thisptr))

    def add_product(self, Species sp, coeff=None):
        """add_product(sp, coeff=None)

        Append a product to the end.

        Args:
            sp (Species): A new product.
            coeff (Integer): A stoichiometry coefficient.

        """
        if coeff is not None:
            self.thisptr.add_product(deref(sp.thisptr), coeff)
        else:
            self.thisptr.add_product(deref(sp.thisptr))

    def set_reactant_coefficient(self, Integer index, Real coeff):
        """set_reactant_coefficient(index, coeff)

        Set a stoichiometry coefficient of a reactant at the given index.

        Args:
            index (Integer): An index pointing the target reactant.
            coeff (Integer): A stoichiometry coefficient.

        """
        self.thisptr.set_reactant_coefficient(index, coeff)

    def set_product_coefficient(self, Integer index, Real coeff):
        """set_product_coefficient(index, coeff)

        Set a stoichiometry coefficient of a product at the given index.

        Args:
            index (Integer): An index pointing the target product.
            coeff (Integer): A stoichiometry coefficient.

        """
        self.thisptr.set_product_coefficient(index, coeff)

    def set_ratelaw(self, ratelaw_obj):
        """set_ratelaw(ratelaw_obj)

        Bind a ratelaw.

        Args:
            ratelaw_obj (ODERatelaw): A ratelaw

        """
        self.ratelaw = ratelaw_obj
        self.thisptr.set_ratelaw(deref((<ODERatelaw>(ratelaw_obj.as_base())).thisptr))

    def set_ratelaw_massaction(self, ODERatelawMassAction ratelaw_obj):
        """set_ratelaw_massaction(ratelaw_obj)

        Bind a mass action ratelaw. This will be deprecated soon.

        Args:
            ratelaw_obj (ODERatelawMassAction): A ratelaw

        """
        self.ratelaw = ratelaw_obj
        self.thisptr.set_ratelaw(deref(ratelaw_obj.thisptr))

    def has_ratelaw(self):
        """Return if a ratelaw is bound or not."""
        return self.thisptr.has_ratelaw()

    def is_massaction(self):
        """Return if a mass action ratelaw is bound or not."""
        return self.thisptr.is_massaction()

    def reactants(self):
        """List all reactants.

        Return:
            list: A list of reactant ``Species``.

        """
        cdef vector[Cpp_Species] cpp_reactants = self.thisptr.reactants()
        retval = []
        cdef vector[Cpp_Species].iterator it = cpp_reactants.begin()
        while it != cpp_reactants.end():
            retval.append(
                    Species_from_Cpp_Species(<Cpp_Species*>address(deref(it))))
            inc(it)
        return retval

    def reactants_coefficients(self):
        """reactants_coefficients() -> [Integer]

        List all coefficients for reactants.

        Return:
            list: A list of reactant coefficients.

        """
        cdef vector[Real] coefficients = self.thisptr.reactants_coefficients()
        retval = []
        cdef vector[Real].iterator it = coefficients.begin()
        while it != coefficients.end():
            retval.append(deref(it))
            inc(it)
        return retval

    def products(self):
        """List all products.

        Return:
            list: A list of product ``Species``.

        """
        cdef vector[Cpp_Species] cpp_products = self.thisptr.products()
        retval = []
        cdef vector[Cpp_Species].iterator it = cpp_products.begin()
        while it != cpp_products.end():
            retval.append(
                    Species_from_Cpp_Species(<Cpp_Species*>address(deref(it))))
            inc(it)
        return retval

    def products_coefficients(self):
        """products_coefficients() -> [Integer]

        List all coefficients for products.

        Return:
            list: A list of product coefficients.

        """
        cdef vector[Real] coefficients = self.thisptr.products_coefficients()
        retval = []
        cdef vector[Real].iterator it = coefficients.begin()
        while it != coefficients.end():
            retval.append( deref(it) )
            inc(it)
        return retval

    def as_string(self):
        """as_string() -> str

        Return an unicode string describing this object.

        Returns:
            str: An unicode string describing this object.

        """
        # reactants = self.reactants()
        # reactants_coeff = self.reactants_coefficients()
        # products = self.products()
        # products_coeff = self.products_coefficients()
        # leftside = ""
        # rightside = ""
        # retval = ""
        # first = True
        # for (sp, coeff) in zip(reactants, reactants_coeff):
        #     s = "{0}({1})".format(coeff, sp.serial())
        #     if first == True:
        #         leftside = s
        #         first = False
        #     else:
        #         leftside = "{} + {}".format(leftside, s)
        # first = True
        # for (sp, coeff) in zip(products, products_coeff):
        #     s = "{0}({1})".format(coeff, sp.serial())
        #     if first == True:
        #         rightside += s
        #         first = False
        #     else:
        #         rightside = "{} + {}".format(retval, s)
        # s = ""
        # if self.has_ratelaw():
        #     s = "HAVE"
        # else:
        #     s = "DON'T HAVE"
        # if self.is_massaction():
        #     k_desc = "k = {:f}\t {} Ratelaw".format(self.k(), s)
        # else:
        #     k_desc = "\t {} Ratelaw".format(s)
        # retval = "{} ---> {}\t{}".format(leftside, rightside, k_desc)
        # return retval
        return self.thisptr.as_string().decode('UTF-8')

cdef ODEReactionRule ODEReactionRule_from_Cpp_ODEReactionRule(Cpp_ODEReactionRule *s):
    cdef Cpp_ODEReactionRule *new_obj = new Cpp_ODEReactionRule(deref(s))
    ret = ODEReactionRule()
    del ret.thisptr
    ret.thisptr = new_obj
    return ret

cdef class ODENetworkModel:
    """A network model class for ODE simulations.

    ODENetworkModel(NetworkModel m=None)

    """

    def __init__(self, NetworkModel m = None):
        """Constructor.

        Args:
            m (NetworkModel, optional): A network model.

        """
        pass

    def __cinit__(self, NetworkModel m = None):
        # self.thisptr = new shared_ptr[Cpp_ODENetworkModel](
        #     <Cpp_ODENetworkModel*>(new Cpp_ODENetworkModel()))
        if m == None:
            self.thisptr = new shared_ptr[Cpp_ODENetworkModel](
                <Cpp_ODENetworkModel*>(new Cpp_ODENetworkModel()))
        else:
            self.thisptr = new shared_ptr[Cpp_ODENetworkModel](
                (<Cpp_ODENetworkModel*>(new Cpp_ODENetworkModel(deref(m.thisptr)))))

    def __dealloc__(self):
        del self.thisptr

    def update_model(self):
        """Update self to fit the given NetworkModel."""
        self.thisptr.get().update_model()

    def has_network_model(self):
        """Return if this model is bound to a NetworkModel or not."""
        return self.thisptr.get().has_network_model()

    def ode_reaction_rules(self):
        """ode_reaction_rules() -> [ODEReactionRule]

        Return a list of ODE reaction rules.

        """
        cdef vector[Cpp_ODEReactionRule] cpp_rules = self.thisptr.get().ode_reaction_rules()
        retval = []
        cdef vector[Cpp_ODEReactionRule].iterator it = cpp_rules.begin()
        while it != cpp_rules.end():
            retval.append(ODEReactionRule_from_Cpp_ODEReactionRule(address(deref(it))))
            inc(it)
        return retval

    def num_reaction_rules(self):
        """Return a number of reaction rules contained in the model."""
        return self.thisptr.get().num_reaction_rules()

    def add_reaction_rule(self, rr):
        """add_reaction_rule(rr)

        Add a new reaction rule.

        Args:
            rr (ReactionRule or ODEReactionRule): A new reaction rule.

        """
        if isinstance(rr, ODEReactionRule):
            self.thisptr.get().add_reaction_rule(deref((<ODEReactionRule>rr).thisptr))
        elif isinstance(rr, ReactionRule):
            self.thisptr.get().add_reaction_rule(deref((<ReactionRule>rr).thisptr))
        else:
            raise ValueError("invalid argument {}".format(repr(rr)))

    def list_species(self):
        """Return a list of species, contained in reaction rules in the model."""
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
    """ A class running the simulation with the ode algorithm.

    ODESimulator(m, w, solver_type)

    """

    def __init__(self, arg1, arg2 = None, arg3 = None):
        """Constructor.

        Args:
            m (ODENetworkModel or NetworkModel): A model
            w (LatticeWorld): A world
            solver_type (int, optional): a type of the ode solver.
                Choose one from RUNGE_KUTTA_CASH_KARP54, ROSENBROCK4 and EULER.

        """
        pass

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
        """Initialize the simulator."""
        self.thisptr.initialize()

    def step(self, upto = None):
        """step(upto=None) -> bool

        Step the simulation.

        Args:
            upto (Real, optional): the time which to step the simulation up to

        Returns:
            bool: True if the simulation did not reach the given time.
                When upto is not given, nothing will be returned.

        """
        if upto is None:
            self.thirptr.step()
        else:
            return self.thisptr.step(upto)

    def next_time(self):
        """Return the scheduled time for the next step."""
        return self.thisptr.next_time()

    def t(self):
        """Return the time."""
        return self.thisptr.t()

    def set_t(self, Real t_new):
        """set_t(t)

        Set the current time.

        Args:
            t (Real): a current time.

        """
        self.thisptr.set_t(t_new)

    def dt(self):
        """Return the step interval."""
        return self.thisptr.dt()

    def set_dt(self, dt_new):
        """set_dt(dt)

        Set a step interval.

        Args:
            dt (Real): a step interval

        """
        self.thisptr.set_dt(dt_new)

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.num_steps()

    def absolute_tolerance(self):
        """Return the absolute tolerance."""
        return self.thisptr.absolute_tolerance()

    def set_absolute_tolerance(self, Real abs_tol):
        """set_absolute_tolerance(abs_tol)

        Set the absolute tolerance.

        Args:
            abs_tol (Real): an absolute tolerance.

        """
        self.thisptr.set_absolute_tolerance(abs_tol)

    def relative_tolerance(self):
        """Return the relative tolerance."""
        return self.thisptr.relative_tolerance()

    def set_relative_tolerance(self, Real rel_tol):
        """set_relative_tolerance(rel_tol)

        Set the relative tolerance.

        Args:
            rel_tol (Real): an relative tolerance.

        """
        self.thisptr.set_relative_tolerance(rel_tol)

    def run(self, Real duration, observers=None):
        """run(duration, observers)

        Run the simulation.

        Args:
            duration (Real): a duration for running a simulation.
                A simulation is expected to be stopped at t() + duration.
            observers (list of Obeservers, optional): observers

        """
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
    """ A factory class creating a ODEWorld instance and a ODESimulator instance.

    ODEFactory(solvertype=None, dt=None)

    """

    def __init__(self, solvertype = None, dt = None):
        """Constructor.

        Args:
            solvertype (int, optional): a type of the ode solver.
                Choose one from RUNGE_KUTTA_CASH_KARP54, ROSENBROCK4 and EULER.
            dt (Real, optional): a default step interval.
        """
        pass

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
        """create_world(arg1=None) -> ODEWorld

        Return a ODEWorld instance.

        Args:
            arg1 (Real3): The lengths of edges of a ODEWorld created

            or

            arg1 (str): The path of a HDF5 file for ODEWorld

        Returns:
            ODEWorld: the created world

        """
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
        """create_simulator(arg1, arg2) -> ODESimulator

        Return a ODESimulator instance.

        Args:
            arg1 (ODEWorld): a world

            or

            arg1 (ODENetworkModel or NetworkModel): a simulation model
            arg2 (ODEWorld): a world

        Returns:
            ODESimulator: the created simulator

        """
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

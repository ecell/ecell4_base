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

    def __init__(self, edge_lengths=None):
        """Constructor.

        Parameters
        ----------
        edge_lengths : Real3, optional
            A size of the World.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, edge_lengths=None):
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

    def actual_lengths(self):
        """actual_lengths() -> Real3

        Return the actual edge lengths of the world.
        Same as ``edge_lengths``.
        """
        cdef Cpp_Real3 lengths = self.thisptr.get().actual_lengths()
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

        Parameters
        ----------
        sp : Species, optional
            a species whose molecules you count

        Returns
        -------
        Integer:
            the number of molecules (of a given species)

        """
        return self.thisptr.get().num_molecules(deref(sp.thisptr))

    def num_molecules_exact(self, Species sp):
        """num_molecules_exact(sp) -> Integer

        Return the number of molecules of a given species.
        A value is rounded to an integer. See get_value_exact also.

        Parameters
        ----------
        sp : Species
            a species whose molecules you count

        Returns
        -------
        Integer:
            the number of molecules of a given species

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

    def new_particle(self, arg1, Real3 arg2=None):
        """new_particle(arg1, arg2=None) -> (ParticleID, Particle)

        Create a new particle.

        Parameters
        ----------
        arg1 : Particle
            A particle to be placed.

        or

        arg1 : Species
            A species of a particle
        arg2 : Real3
            A position to place a particle

        Returns
        -------
        tuple:
            A pair of ParticleID and Particle of a new particle

        """
        cdef pair[pair[Cpp_ParticleID, Cpp_Particle], bool] retval

        if arg2 is None:
            retval = self.thisptr.get().new_particle(deref((<Particle> arg1).thisptr))
        else:
            retval = self.thisptr.get().new_particle(deref((<Species> arg1).thisptr), deref(arg2.thisptr))
        return ((ParticleID_from_Cpp_ParticleID(address(retval.first.first)), Particle_from_Cpp_Particle(address(retval.first.second))), retval.second)

    def add_molecules(self, Species sp, Integer num, shape=None):
        """add_molecules(sp, num, shape=None)

        Add some molecules.

        Parameters
        ----------
        sp : Species
            a species of molecules to add
        num : Integer
            the number of molecules to add
        shape : Shape, optional
            a shape to add molecules on [not supported yet]

        """
        if shape is None:
            self.thisptr.get().add_molecules(deref(sp.thisptr), num)
        else:
            self.thisptr.get().add_molecules(
                deref(sp.thisptr), num, deref((<Shape>(shape.as_base())).thisptr))

    def remove_molecules(self, Species sp, Integer num):
        """remove_molecules(sp, num)

        Remove molecules

        Parameters
        ----------
        sp : Species
            a species whose molecules to remove
        num : Integer
            a number of molecules to be removed

        """
        self.thisptr.get().remove_molecules(deref(sp.thisptr), num)

    def get_value(self, Species sp):
        """get_value(sp) -> Real

        Return the value matched to a given species.

        Parameters
        ----------
        sp : Species
            a pattern whose value you get

        Returns
        -------
        Real:
            the value matched to a given species

        """
        return self.thisptr.get().get_value(deref(sp.thisptr))

    def get_value_exact(self, Species sp):
        """get_value_exact(sp) -> Real

        Return the value connected to a given species.

        Parameters
        ----------
        sp : Species
            a species whose value you get

        Returns
        -------
        Real:
            the value connected to a given species

        """
        return self.thisptr.get().get_value(deref(sp.thisptr))

    def set_value(self, Species sp, Real value):
        """set_value(sp, value)

        Set the value of the given species.

        Parameters
        ----------
        sp : Species
            a species whose value you set
        value : Real
            a value set

        """
        self.thisptr.get().set_value(deref(sp.thisptr), value)

    def save(self, filename):
        """save(filename)

        Save the current state to a HDF5 file.

        Parameters
        ----------
        filename : str
            a file name to be saved.

        """
        self.thisptr.get().save(tostring(filename))

    def load(self, filename):
        """load(filename)

        Load a HDF5 file to the current state.

        Parameters
        ----------
        filename : str
            a file name to be loaded.

        """
        self.thisptr.get().load(tostring(filename))

    def has_species(self, Species sp):
        """has_species(sp) -> bool

        Check if the given species is belonging to this.

        Parameters
        ----------
        sp : Species
            a species to be checked.

        Returns
        -------
        bool:
            True if the given species is contained.

        """
        return self.thisptr.get().has_species(deref(sp.thisptr))

    def reserve_species(self, Species sp):
        """reserve_species(sp)

        Reserve a value for the given species. Use set_value.

        Parameters
        ----------
        sp : Species
            a species to be reserved.

        """
        self.thisptr.get().reserve_species(deref(sp.thisptr))

    def release_species(self, Species sp):
        """release_species(sp)

        Release a value for the given species.
        This function is mainly for developers.

        Parameters
        ----------
        sp : Species
            a species to be released.

        """
        self.thisptr.get().release_species(deref(sp.thisptr))

    def bind_to(self, m):
        """bind_to(m)

        Bind a model.

        Parameters
        ----------
        m : ODENetworkModel or NetworkModel
            a model to be bound

        """
        if isinstance(m, ODENetworkModel):
            self.thisptr.get().bind_to(deref((<ODENetworkModel>m).thisptr))
        else:
            self.thisptr.get().bind_to(Cpp_Model_from_Model(m))

    def evaluate(self, rr):
        if isinstance(rr, ReactionRule):
            return self.thisptr.get().evaluate(deref((<ReactionRule>rr).thisptr))
        elif isinstance(rr, ODEReactionRule):
            return self.thisptr.get().evaluate(deref((<ODEReactionRule>rr).thisptr))
        else:
            raise ValueError(
                "A ReactionRule or ODEReactionRule must be given [{}].".format(repr(rr)))

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

    def as_string(self):
        """"Return a name of the function"""
        return self.thisptr.get().as_string().decode('UTF-8')

    def as_base(self):
        """Return self as a base class. Only for developmental use."""
        return self

    def to_derivative(self):
        r = ODERatelawMassAction_from_Cpp_ODERatelaw(deref(self.thisptr) )
        if r !=  None:
            return r

        r = ODERatelawCallback_from_Cpp_ODERatelaw(deref(self.thisptr) )
        if r != None:
            return r

        raise ValueError("Invalid Ratelaw Type")
    def __reduce__(self):
        return self.to_derivative().__reduce__()


cdef ODERatelaw ODERatelaw_from_Cpp_ODERatelaw(shared_ptr[Cpp_ODERatelaw] s):
    r = ODERatelaw()
    r.thisptr.swap(s)
    return r

cdef ODERatelawMassAction ODERatelawMassAction_from_Cpp_ODERatelaw(shared_ptr[Cpp_ODERatelaw] s):
    r = ODERatelawMassAction(0.01)
    cdef shared_ptr[Cpp_ODERatelawMassAction] temp = to_ODERatelawMassAction(s)
    if temp.get() == NULL:
        return None
    r.thisptr.swap(temp)
    return r

cdef ODERatelawCallback ODERatelawCallback_from_Cpp_ODERatelaw(shared_ptr[Cpp_ODERatelaw] s):
    r = ODERatelawCallback(lambda x:x)
    cdef shared_ptr[Cpp_ODERatelawCythonCallback] temp = to_ODERatelawCythonCallback(s)
    if temp.get() == NULL:
        return None
    r.thisptr.swap(temp)
    return r

cdef class ODERatelawMassAction:
    """A class for mass action ratelaws.

    ODERatelawMassAction(Real k)

    """

    def __init__(self, Real k):
        """Constructor.

        Parameters
        ----------
        k : Real
            A kinetic rate for the mass action.

        """
        pass

    def __cinit__(self, Real k):
        self.thisptr = new shared_ptr[Cpp_ODERatelawMassAction](
                <Cpp_ODERatelawMassAction*>(new Cpp_ODERatelawMassAction(k)))

    def __dealloc__(self):
        del self.thisptr

    def is_available(self):
        """Check if this ratelaw is available or not. Return True always."""
        return self.thisptr.get().is_available()

    def set_k(self, Real k):
        """set_k(k)

        Set a kinetic rate constant.

        Parameters
        ----------
        k : float
            A kinetic rate constant.

        """
        #self.get().thisptr.set_k(k)
        self.thisptr.get().set_k(k)

    def get_k(self):
        """Return the kinetic rate constant as a float value."""
        #return self.get().thisptr.get_k()
        return self.thisptr.get().get_k()

    def as_string(self):
        """"Return a name of the function"""
        return self.thisptr.get().as_string().decode('UTF-8')

    def as_base(self):
        """Return self as a base class. Only for developmental use."""
        base_type = ODERatelaw()
        del base_type.thisptr
        base_type.thisptr = new shared_ptr[Cpp_ODERatelaw](
                <shared_ptr[Cpp_ODERatelaw]>(deref(self.thisptr)))
        return base_type
    def __reduce__(self):
        return (__rebuild_ode_ratelaw, ("ODERatelawMassAction", self.as_string(), self.get_k() ) )


# cdef double indirect_function(
cdef indirect_function(
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

    ODERatelawCallback(pyfunc, name)

    """

    def __init__(self, pyfunc, name=None):
        """Constructor.

        Parameters
        ----------
        pyfunc : function
            A Python function for the callback.
            See set_callback function of this class for details.
        name : string, optional
            A name of the function

        """
        pass

    def __cinit__(self, pyfunc, name=None):
        if name is None:
            self.thisptr = new shared_ptr[Cpp_ODERatelawCythonCallback](
                <Cpp_ODERatelawCythonCallback*>(new Cpp_ODERatelawCythonCallback(
                    <Stepladder_Functype>indirect_function, <void*>pyfunc,
                    <OperateRef_Functype>inc_ref, <OperateRef_Functype>dec_ref,
                    tostring(pyfunc.__name__))))
        else:
            self.thisptr = new shared_ptr[Cpp_ODERatelawCythonCallback](
                <Cpp_ODERatelawCythonCallback*>(new Cpp_ODERatelawCythonCallback(
                    <Stepladder_Functype>indirect_function, <void*>pyfunc,
                    <OperateRef_Functype>inc_ref, <OperateRef_Functype>dec_ref,
                    tostring(name))))
        self.pyfunc = pyfunc

    def __dealloc__(self):
        del self.thisptr

    def set_callback(self, pyfunc):
        """set_callback(pyfunc)

        Parameters
        ----------
        pyfunc : function
            A Python function for the callback
            The function must accept five arguments, and return a velocity.
            The number of reactants, the number of products, a volume,
            the current time, and a ODEReactionRule are given as the
            arguments in this order.

        Examples
        --------
        The following callback represents a simple Michaelis-Menten-like
        equation:

        >>> rl = ODERatelawCallback()
        >>> rl.set_callback(lambda r, p, v, t, rr: 2.0 * r[0] * r[1] / (1.0 + r[1]))

        Here, we expect that the first reactant is an enzyme,
        and that the second one is a substrate.

        """
        self.thisptr.get().set_callback_pyfunc(<Python_CallbackFunctype>pyfunc)
        self.pyfunc = pyfunc
    def get_callback(self):
        return <object>self.thisptr.get().get_callback_pyfunc()

    def set_name(self, name):
        """"Set the name of a function"""
        self.thisptr.get().set_name(tostring(name))

    def as_string(self):
        """"Return a name of the function"""
        return self.thisptr.get().as_string().decode('UTF-8')

    def as_base(self):
        """Return self as a base class. Only for developmental use."""
        retval = ODERatelaw()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_ODERatelaw](
            <shared_ptr[Cpp_ODERatelaw]>deref(self.thisptr))
        return retval
    def get_pyfunc(self):
        return self.pyfunc

    def __reduce__(self):
        import sys
        loaded_modules = sys.modules.keys()
        if not  "dill" in loaded_modules:
            raise RuntimeError("dill module is required for pickling user-defined function")
        return (__rebuild_ode_ratelaw, ("ODERatelawCallback", self.as_string(), self.get_callback()) )

def __rebuild_ode_ratelaw(ratelaw_type, name, param):
    if ratelaw_type == "ODERatelawCallback":
        m = ODERatelawCallback(param, name)
        return m
    elif ratelaw_type == "ODERatelawMassAction":
        m = ODERatelawMassAction(param)
        #m.set_name(name)
        return m
    else:
        raise ValueError("Invalid Ratelaw Type")
    

cdef class ODEReactionRule:
    """A class representing a reaction rule between ``Species``, which accepts at most
    one rate law to calculate the flux.

    ODEReactionRule(rr)

    """

    def __init__(self, *args):
        """Constructor.

        Parameters
        ----------
        rr : ReactionRule

        """
        pass

    def __cinit__(self, *args):
        if len(args) == 0:
            self.thisptr = new Cpp_ODEReactionRule()
            self.ratelaw = None
        elif len(args) == 1 and isinstance(args[0], ReactionRule):
            self.thisptr = new Cpp_ODEReactionRule(deref((<ReactionRule>args[0]).thisptr))
            self.ratelaw = None
        else:
            raise ValueError("The invalid arguments are given.")

    def __dealloc__(self):
        del self.thisptr

    def k(self):
        """Return the kinetic rate constant as a float value."""
        return self.thisptr.k()

    def set_k(self, Real k):
        """set_k(k)

        Set a kinetic rate constant.

        Parameters
        ----------
        k : float
            A kinetic rate constant.

        """
        self.thisptr.set_k(k)

    def add_reactant(self, Species sp, coeff=None):
        """add_reactant(sp, coeff=None)

        Append a reactant to the end.

        Parameters
        ----------
        sp : Species
            A new reactant.
        coeff : Integer
            A stoichiometry coefficient.

        """
        if coeff is not None:
            self.thisptr.add_reactant(deref(sp.thisptr), coeff)
        else:
            self.thisptr.add_reactant(deref(sp.thisptr))

    def add_product(self, Species sp, coeff=None):
        """add_product(sp, coeff=None)

        Append a product to the end.

        Parameters
        ----------
        sp : Species
            A new product.
        coeff : Integer
            A stoichiometry coefficient.

        """
        if coeff is not None:
            self.thisptr.add_product(deref(sp.thisptr), coeff)
        else:
            self.thisptr.add_product(deref(sp.thisptr))

    def set_reactant_coefficient(self, Integer index, Real coeff):
        """set_reactant_coefficient(index, coeff)

        Set a stoichiometry coefficient of a reactant at the given index.

        Parameters
        ----------
        index : Integer
            An index pointing the target reactant.
        coeff : Integer
            A stoichiometry coefficient.

        """
        self.thisptr.set_reactant_coefficient(index, coeff)

    def set_product_coefficient(self, Integer index, Real coeff):
        """set_product_coefficient(index, coeff)

        Set a stoichiometry coefficient of a product at the given index.

        Parameters
        ----------
        index : Integer
            An index pointing the target product.
        coeff : Integer
            A stoichiometry coefficient.

        """
        self.thisptr.set_product_coefficient(index, coeff)

    def set_ratelaw(self, ratelaw_obj):
        """set_ratelaw(ratelaw_obj)

        Bind a ratelaw.

        Parameters
        ----------
        ratelaw_obj : ODERatelaw
            A ratelaw

        """
        self.ratelaw = ratelaw_obj
        self.thisptr.set_ratelaw(deref((<ODERatelaw>(ratelaw_obj.as_base())).thisptr))

    def set_ratelaw_massaction(self, ODERatelawMassAction ratelaw_obj):
        """set_ratelaw_massaction(ratelaw_obj)

        Bind a mass action ratelaw. This will be deprecated soon.

        Parameters
        ----------
        ratelaw_obj : ODERatelawMassAction
            A ratelaw

        """
        self.ratelaw = ratelaw_obj
        self.thisptr.set_ratelaw(deref(ratelaw_obj.thisptr))

    def has_ratelaw(self):
        """Return if a ratelaw is bound or not."""
        return self.thisptr.has_ratelaw()

    def get_ratelaw(self):
        """Return a ratelaw"""
        return ODERatelaw_from_Cpp_ODERatelaw(self.thisptr.get_ratelaw())

    def is_massaction(self):
        """Return if a mass action ratelaw is bound or not."""
        return self.thisptr.is_massaction()

    def reactants(self):
        """List all reactants.

        Returns
        -------
        list:
            A list of reactant ``Species``.

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

        Returns
        -------
        list:
            A list of reactant coefficients.

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

        Returns
        -------
        list:
            A list of product ``Species``.

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

        Returns
        -------
        list:
            A list of product coefficients.

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

        Returns
        -------
        str:
            An unicode string describing this object.

        """
        return self.thisptr.as_string().decode('UTF-8')

    def __reduce__(self):
        if self.has_ratelaw():
            ratelaw = self.get_ratelaw()
        else:
            ratelaw = None
        return (__rebuild_ode_reaction_rule, (self.reactants(), self.products(), self.reactants_coefficients(), self.products_coefficients(), ratelaw))

def __rebuild_ode_reaction_rule(reactants, products, reactants_coefficients, products_coefficients, ratelaw):
    rr = ODEReactionRule()
    for sp, coef in zip(reactants, reactants_coefficients):
        rr.add_reactant(sp, coef)
    for sp, coef in zip(products, products_coefficients):
        rr.add_product(sp, coef)

    if ratelaw is None:
        pass
    else:
        rr.set_ratelaw(ratelaw)
    return rr

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

    def __init__(self, m=None):
        """Constructor.

        Parameters
        ----------
        m : Model, optional
            A network model.

        """
        pass

    def __cinit__(self, m=None):
        # self.thisptr = new shared_ptr[Cpp_ODENetworkModel](
        #     <Cpp_ODENetworkModel*>(new Cpp_ODENetworkModel()))
        if m == None:
            self.thisptr = new shared_ptr[Cpp_ODENetworkModel](
                <Cpp_ODENetworkModel*>(new Cpp_ODENetworkModel()))
        # else:
        #     # self.thisptr = new shared_ptr[Cpp_ODENetworkModel](
        #     #     (<Cpp_ODENetworkModel*>(new Cpp_ODENetworkModel(deref(m.thisptr)))))
        #     self.thisptr = new shared_ptr[Cpp_ODENetworkModel](
        #         (<Cpp_ODENetworkModel*>(new Cpp_ODENetworkModel(m.thisptr))))
        elif isinstance(m, Model):
            self.thisptr = new shared_ptr[Cpp_ODENetworkModel](
                (<Cpp_ODENetworkModel*>(new Cpp_ODENetworkModel((<Model>m).thisptr))))
        elif isinstance(m, NetworkModel):
            self.thisptr = new shared_ptr[Cpp_ODENetworkModel](
                (<Cpp_ODENetworkModel*>(new Cpp_ODENetworkModel(<shared_ptr[Cpp_Model]>((<NetworkModel>m).thisptr)))))
        elif isinstance(m, NetfreeModel):
            self.thisptr = new shared_ptr[Cpp_ODENetworkModel](
                (<Cpp_ODENetworkModel*>(new Cpp_ODENetworkModel(<shared_ptr[Cpp_Model]>((<NetfreeModel>m).thisptr)))))
        else:
            raise ValueError('Unsupported model type was given.')

    def __dealloc__(self):
        del self.thisptr

    def update_model(self):
        """Update self to fit the given Model."""
        self.thisptr.get().update_model()

    def has_network_model(self):
        """Return if this model is bound to a Model or not."""
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

    def reaction_rules(self):
        """reaction_rules() -> [ODEReactionRule]

        Return a list of ODE reaction rules.

        """
        cdef vector[Cpp_ODEReactionRule] cpp_rules = self.thisptr.get().reaction_rules()
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

        Parameters
        ----------
        rr : ReactionRule or ODEReactionRule
            A new reaction rule.

        """
        if isinstance(rr, ODEReactionRule):
            self.thisptr.get().add_reaction_rule(deref((<ODEReactionRule>rr).thisptr))
        elif isinstance(rr, ReactionRule):
            self.thisptr.get().add_reaction_rule(deref((<ReactionRule>rr).thisptr))
        else:
            raise ValueError("invalid argument {}".format(repr(rr)))

    def add_reaction_rules(self, rrs):
        if isinstance(rrs, list):
            for rr in rrs:
                self.add_reaction_rule(rr)
        else:
            self.add_reaction_rule(rrs)

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

    def __reduce__(self):
        return (__rebuild_ode_network_model, (self.reaction_rules(), ))

def __rebuild_ode_network_model(rrs):
    m = ODENetworkModel()
    m.add_reaction_rules(rrs)
    return m

cdef ODENetworkModel ODENetworkModel_from_Cpp_ODENetworkModel(
    shared_ptr[Cpp_ODENetworkModel] m):
    r = ODENetworkModel()
    r.thisptr.swap(m)
    return r

# ODESolverType:
(
    RUNGE_KUTTA_CASH_KARP54,
    ROSENBROCK4_CONTROLLER,
    EULER,
) = (0, 1, 2)

cdef Cpp_ODESolverType translate_solver_type(solvertype_constant):
    if solvertype_constant == RUNGE_KUTTA_CASH_KARP54:
        return RUNGE_KUTTA_CASH_KARP54
    elif solvertype_constant == ROSENBROCK4_CONTROLLER:
        return Cpp_ROSENBROCK4_CONTROLLER
    elif solvertype_constant == EULER:
        return Cpp_EULER
    else:
        raise ValueError(
            "invalid solver type was given [{0}]".format(repr(solvertype_constant)))

# ODERatelawType:
(
    ABSTRACT_TYPE,
    MASSACTION_TYPE,
    PYTHON_CALLBACK_TYPE,
    CPP_CALLBACK_TYPE,
) = (0, 1, 2, 3)


cdef class ODESimulator:
    """ A class running the simulation with the ode algorithm.

    ODESimulator(m, w, solver_type)

    """

    def __init__(self, arg1, arg2=None, arg3=None):
        """Constructor.

        Parameters
        ----------
        m : ODENetworkModel or Model
            A model
        w : ODEWorld
            A world
        solver_type : int, optional
            a type of the ode solver.
            Choose one from RUNGE_KUTTA_CASH_KARP54, ROSENBROCK4_CONTROLLER and EULER.

        """
        pass

    def __cinit__(self, arg1, arg2=None, arg3=None):
        if arg2 is None or not isinstance(arg2, ODEWorld):
            if not isinstance(arg1, ODEWorld):
                raise ValueError(
                    "An invalid value [{}] for the first argument.".format(repr(arg1))
                    + " ODEWorld is needed.")

            if arg2 is None:
                self.thisptr = new Cpp_ODESimulator(
                    deref((<ODEWorld>arg1).thisptr))
            else:
                self.thisptr = new Cpp_ODESimulator(
                    deref((<ODEWorld>arg1).thisptr),
                    translate_solver_type(arg2))
        else:
            if isinstance(arg1, ODENetworkModel):
                if arg3 is None:
                    self.thisptr = new Cpp_ODESimulator(
                        deref((<ODENetworkModel>arg1).thisptr),
                        deref((<ODEWorld>arg2).thisptr))
                else:
                    self.thisptr = new Cpp_ODESimulator(
                        deref((<ODENetworkModel>arg1).thisptr),
                        deref((<ODEWorld>arg2).thisptr),
                        translate_solver_type(arg3))
            # elif isinstance(arg1, Model):
            else:
                if arg3 is None:
                    self.thisptr = new Cpp_ODESimulator(
                        Cpp_Model_from_Model(arg1),
                        deref((<ODEWorld>arg2).thisptr))
                else:
                    self.thisptr = new Cpp_ODESimulator(
                        Cpp_Model_from_Model(arg1),
                        deref((<ODEWorld>arg2).thisptr),
                        translate_solver_type(arg3))
            # else:
            #     raise ValueError(
            #         "An invalid value [{}] for the first argument.".format(repr(arg1))
            #         + " ODENetworkModel or Model is needed.")

    def __dealloc__(self):
        del self.thisptr

    def initialize(self):
        """Initialize the simulator."""
        self.thisptr.initialize()

    def step(self, upto=None):
        """step(upto=None) -> bool

        Step the simulation.

        Parameters
        ----------
        upto : Real, optional
            the time which to step the simulation up to

        Returns
        -------
        bool:
            True if the simulation did not reach the given time.
            When upto is not given, nothing will be returned.

        """
        if upto is None:
            self.thisptr.step()
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

        Parameters
        ----------
        t : Real
            a current time.

        """
        self.thisptr.set_t(t_new)

    def dt(self):
        """Return the step interval."""
        return self.thisptr.dt()

    def set_dt(self, dt_new):
        """set_dt(dt)

        Set a step interval.

        Parameters
        ----------
        dt : Real
            a step interval

        """
        self.thisptr.set_dt(dt_new)

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.num_steps()

    def check_reaction(self):
        """Return if any reaction occurred at the last step, or not.
        This function always returns False."""
        return self.thisptr.check_reaction()

    def absolute_tolerance(self):
        """Return the absolute tolerance."""
        return self.thisptr.absolute_tolerance()

    def set_absolute_tolerance(self, Real abs_tol):
        """set_absolute_tolerance(abs_tol)

        Set the absolute tolerance.

        Parameters
        ----------
        abs_tol : Real
            an absolute tolerance.

        """
        self.thisptr.set_absolute_tolerance(abs_tol)

    def relative_tolerance(self):
        """Return the relative tolerance."""
        return self.thisptr.relative_tolerance()

    def set_relative_tolerance(self, Real rel_tol):
        """set_relative_tolerance(rel_tol)

        Set the relative tolerance.

        Parameters
        ----------
        rel_tol : Real
            an relative tolerance.

        """
        self.thisptr.set_relative_tolerance(rel_tol)

    def model(self):
        """Return the model bound."""
        return ODENetworkModel_from_Cpp_ODENetworkModel(self.thisptr.model())

    def world(self):
        """Return the world bound."""
        return ODEWorld_from_Cpp_ODEWorld(self.thisptr.world())

    def run(self, Real duration, observers=None):
        """run(duration, observers)

        Run the simulation.

        Parameters
        ----------
        duration : Real
            a duration for running a simulation.
                A simulation is expected to be stopped at t() + duration.
        observers : list of Obeservers, optional
            observers

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

    ODEFactory(ODESolverType solver_type=None, Real dt=None, Real abs_tol=None, Real rel_tol=None)

    """

    def __init__(self, solver_type=None, dt=None, abs_tol=None, rel_tol=None):
        """Constructor.

        Parameters
        ----------
        solver_type : int, optional
            a type of the ode solver.
            Choose one from RUNGE_KUTTA_CASH_KARP54, ROSENBROCK4_CONTROLLER and EULER.
        dt : Real, optional
            a default step interval.
        abs_tol : Real, optional
            absolute tolerance.
        rel_tol : Real, optional
            relative tolerance.

        """
        pass

    def __cinit__(self, solver_type=None, dt=None, abs_tol=None, rel_tol=None):
        self.thisptr = new Cpp_ODEFactory(
            Cpp_ODEFactory.default_solver_type() if solver_type is None else translate_solver_type(solver_type),
            Cpp_ODEFactory.default_dt() if dt is None else <Real>dt,
            Cpp_ODEFactory.default_abs_tol() if abs_tol is None else <Real>abs_tol,
            Cpp_ODEFactory.default_rel_tol() if rel_tol is None else <Real>rel_tol)

    def rng(self, GSLRandomNumberGenerator rng):
        """rng(GSLRandomNumberGenerator) -> ODEFactory

        Just return self. This method is for the compatibility between Factory classes.

        """
        cdef Cpp_ODEFactory *ptr = self.thisptr.rng_ptr(deref(rng.thisptr))
        assert ptr == self.thisptr
        return self

    def __dealloc__(self):
        del self.thisptr

    def create_world(self, arg1=None):
        """create_world(arg1=None) -> ODEWorld

        Return a ODEWorld instance.

        Parameters
        ----------
        arg1 : Real3
            The lengths of edges of a ODEWorld created

        or

        arg1 : str
            The path of a HDF5 file for ODEWorld

        Returns
        -------
        ODEWorld:
            the created world

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
                shared_ptr[Cpp_ODEWorld](self.thisptr.create_world(tostring(arg1))))
        raise ValueError("invalid argument")

    # def create_simulator(self, arg1, ODEWorld arg2=None):
    #     if arg2 is None:
    #         return ODESimulator_from_Cpp_ODESimulator(
    #             self.thisptr.create_simulator(deref((<ODEWorld>arg1).thisptr)))
    #     else:
    #         return ODESimulator_from_Cpp_ODESimulator(
    #             self.thisptr.create_simulator(
    #                 deref((<ODENetworkModel>arg1).thisptr), deref(arg2.thisptr)))

    def create_simulator(self, arg1, arg2=None):
        """create_simulator(arg1, arg2) -> ODESimulator

        Return a ODESimulator instance.

        Parameters
        ----------
        arg1 : ODEWorld
            a world

        or

        arg1 : ODENetworkModel or NetworkModel
            a simulation model
        arg2 : ODEWorld
            a world

        Returns
        -------
        ODESimulator:
            the created simulator

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
            else: # elif isinstance(arg1, NetworkModel):
                return ODESimulator_from_Cpp_ODESimulator(
                    self.thisptr.create_simulator(
                        Cpp_Model_from_Model(arg1), # (<NetworkModel>arg1).thisptr,
                        deref((<ODEWorld>arg2).thisptr)))

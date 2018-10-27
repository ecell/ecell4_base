import collections
from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address
from libcpp.string cimport string
from libcpp.vector cimport vector

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *

from cpython cimport PyObject, Py_XINCREF, Py_XDECREF

from deprecated import deprecated


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

    @deprecated(suggest="edge_lengths()")
    def actual_lengths(self):
        return self.edge_lengths()

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
        m : Model
            a model to be bound

        """
        self.thisptr.get().bind_to(Cpp_Model_from_Model(m))

    def evaluate(self, rr):
        if isinstance(rr, ReactionRule):
            return self.thisptr.get().evaluate(deref((<ReactionRule>rr).thisptr))
        else:
            raise ValueError(
                "A ReactionRule must be given [{}].".format(repr(rr)))

    def as_base(self):
        """Return self as a base class. Only for developmental use."""
        retval = WorldInterface()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_WorldInterface](
            <shared_ptr[Cpp_WorldInterface]>deref(self.thisptr))
        return retval

cdef ODEWorld ODEWorld_from_Cpp_ODEWorld(
    shared_ptr[Cpp_ODEWorld] w):
    r = ODEWorld(Real3(1, 1, 1))
    r.thisptr.swap(w)
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

cdef class ODESimulator:
    """ A class running the simulation with the ode algorithm.

    ODESimulator(m, w, solver_type)

    """

    def __init__(self, ODEWorld w, m=None, solver_type=None):
        """Constructor.

        Parameters
        ----------
        w : ODEWorld
            A world
        m : Model, optional
            A model
        solver_type : int, optional
            a type of the ode solver.
            Choose one from RUNGE_KUTTA_CASH_KARP54, ROSENBROCK4_CONTROLLER and EULER.

        """
        pass

    def __cinit__(self, ODEWorld w, m=None, solver_type=None):
        if m is None and solver_type is None:
            self.thisptr = new Cpp_ODESimulator(deref(w.thisptr))
        elif m is not None and solver_type is None:
            self.thisptr = new Cpp_ODESimulator(deref(w.thisptr), Cpp_Model_from_Model(m))
        elif m is None and solver_type is not None:
            self.thisptr = new Cpp_ODESimulator(deref(w.thisptr), translate_solver_type(solver_type))
        else:
            self.thisptr = new Cpp_ODESimulator(
                deref(w.thisptr), Cpp_Model_from_Model(m), translate_solver_type(solver_type))

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
        return Model_from_Cpp_Model(self.thisptr.model())

    def world(self):
        """Return the world bound."""
        return ODEWorld_from_Cpp_ODEWorld(self.thisptr.world())

    def run(self, Real duration, observers=None, is_dirty=None):
        """run(duration, observers, is_dirty)

        Run the simulation.

        Parameters
        ----------
        duration : Real
            A duration for running a simulation.
            A simulation is expected to be stopped at t() + duration.
        observers : list of Obeservers, optional
            observers
        is_dirty : bool, default True
            If True, call initialize before running.

        """
        cdef vector[shared_ptr[Cpp_Observer]] tmp

        if is_dirty is None:
            if observers is None:
                self.thisptr.run(duration)
            elif isinstance(observers, collections.Iterable):
                for obs in observers:
                    tmp.push_back(deref((<Observer>(obs.as_base())).thisptr))
                self.thisptr.run(duration, tmp)
            else:
                self.thisptr.run(duration,
                    deref((<Observer>(observers.as_base())).thisptr))
        else:
            if observers is None:
                self.thisptr.run(duration, is_dirty)
            elif isinstance(observers, collections.Iterable):
                for obs in observers:
                    tmp.push_back(deref((<Observer>(obs.as_base())).thisptr))
                self.thisptr.run(duration, tmp, is_dirty)
            else:
                self.thisptr.run(duration,
                    deref((<Observer>(observers.as_base())).thisptr), is_dirty)


cdef ODESimulator ODESimulator_from_Cpp_ODESimulator(Cpp_ODESimulator* s):
    r = ODESimulator(
        Model_from_Cpp_Model(s.model()),
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

    def world(self, arg1=None):
        """world(arg1=None) -> ODEWorld

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
                shared_ptr[Cpp_ODEWorld](self.thisptr.world()))
        elif isinstance(arg1, Real3):
            return ODEWorld_from_Cpp_ODEWorld(
                shared_ptr[Cpp_ODEWorld](
                    self.thisptr.world(deref((<Real3>arg1).thisptr))))
        elif isinstance(arg1, str):
            return ODEWorld_from_Cpp_ODEWorld(
                shared_ptr[Cpp_ODEWorld](self.thisptr.world(tostring(arg1))))
        raise ValueError("invalid argument")

    def simulator(self, ODEWorld arg1, arg2=None):
        """simulator(arg1, arg2) -> ODESimulator

        Return a ODESimulator instance.

        Parameters
        ----------
        arg1 : ODEWorld
            a world
        arg2 : Model, optional
            a simulation model

        Returns
        -------
        ODESimulator:
            the created simulator

        """
        if arg2 is None:
            return ODESimulator_from_Cpp_ODESimulator(
                self.thisptr.simulator(deref(arg1.thisptr)))
        else:
            return ODESimulator_from_Cpp_ODESimulator(
                self.thisptr.simulator(
                    deref(arg1.thisptr),
                    Cpp_Model_from_Model(arg2)))

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
        return self.world(arg1)

    def create_simulator(self, ODEWorld arg1, arg2=None):
        """create_simulator(arg1, arg2) -> ODESimulator

        Return a ODESimulator instance.

        Parameters
        ----------
        arg1 : ODEWorld
            a world
        arg2 : Model, optional
            a simulation model

        Returns
        -------
        ODESimulator:
            the created simulator

        """
        return self.simulator(arg1, arg2)

import collections
from cython cimport address
from cython.operator cimport dereference as deref, preincrement as inc
from ecell4.core cimport *
import numbers


## ReactionInfo
cdef class ReactionInfo:
    """A class stores detailed information about a reaction in egfrd.

    ReactionInfo(t, reactants, products)

    """

    def __init__(self, Real t, reactants, products):
        """Constructor.

        Args:
          t (Real): A time when a reaction occurs
          reactants (list): A list of reactants.
            Reactants are given as a pair of ``ParticleID`` and ``Particle``.
          products (list): A list of products.
            Products are given as a pair of ``ParticleID`` and ``Particle``.

        """
        pass  #XXX: only used for doc string

    def __cinit__(self, Real t, reactants, products):
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]] reactants_
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]] products_

        for pid, p in reactants:
            reactants_.push_back(
                pair[Cpp_ParticleID, Cpp_Particle](
                    deref((<ParticleID>pid).thisptr), deref((<Particle>p).thisptr)))
        for pid, p in products:
            products_.push_back(
                pair[Cpp_ParticleID, Cpp_Particle](
                    deref((<ParticleID>pid).thisptr), deref((<Particle>p).thisptr)))

        self.thisptr = new Cpp_ReactionInfo(t, reactants_, products_)

    def __dealloc__(self):
        del self.thisptr

    def t(self):
        """Return a time when a reaction occurred."""
        return self.thisptr.t()

    def reactants(self):
        """Return a list of reactants

        Returns:
            list: A list of pairs of ``ParticleID`` and ``Particle``.

        """
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]] particles
        particles = self.thisptr.reactants()

        retval = []
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]].iterator \
            it = particles.begin()
        while it != particles.end():
            retval.append(
                (ParticleID_from_Cpp_ParticleID(
                     <Cpp_ParticleID*>(address(deref(it).first))),
                 Particle_from_Cpp_Particle(
                     <Cpp_Particle*>(address(deref(it).second)))))
            inc(it)
        return retval

    def products(self):
        """Return a list of products

        Returns:
            list: A list of pairs of ``ParticleID`` and ``Particle``.

        """
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]] particles
        particles = self.thisptr.products()

        retval = []
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]].iterator \
            it = particles.begin()
        while it != particles.end():
            retval.append(
                (ParticleID_from_Cpp_ParticleID(
                     <Cpp_ParticleID*>(address(deref(it).first))),
                 Particle_from_Cpp_Particle(
                     <Cpp_Particle*>(address(deref(it).second)))))
            inc(it)
        return retval

    def __reduce__(self):
        return (ReactionInfo, (self.t(), self.reactants(), self.products()))

cdef ReactionInfo ReactionInfo_from_Cpp_ReactionInfo(Cpp_ReactionInfo* ri):
    cdef Cpp_ReactionInfo *new_obj = new Cpp_ReactionInfo(<Cpp_ReactionInfo> deref(ri))
    r = ReactionInfo(0, [], [])
    del r.thisptr
    r.thisptr = new_obj
    return r

## EGFRDWorld
#  a python wrapper for Cpp_EGFRDWorld
cdef class EGFRDWorld:
    """A class containing the properties of the egfrd world.

    EGFRDWorld(edge_lengths=None, matrix_sizes=None, GSLRandomNumberGenerator rng=None)

    """

    def __init__(self, edge_lengths = None, Integer3 matrix_sizes = None,
                 GSLRandomNumberGenerator rng = None):
        """Constructor.

        Parameters
        ----------
        edge_lengths : Real3, optional
            A size of the World.
        matrix_sizes : Integer3, optional
            A size of a cell matrix.
            The number of cells must be larger than 3, in principle.
        rng : GSLRandomNumberGenerator, optional
            A random number generator.

        """
        pass

    def __cinit__(self, edge_lengths=None, Integer3 matrix_sizes=None,
                  GSLRandomNumberGenerator rng=None):
        cdef string filename

        if rng is not None:
            self.thisptr = new shared_ptr[Cpp_EGFRDWorld](
                new Cpp_EGFRDWorld(
                    deref((<Real3>edge_lengths).thisptr),
                    deref(matrix_sizes.thisptr), deref(rng.thisptr)))
        elif matrix_sizes is not None:
            self.thisptr = new shared_ptr[Cpp_EGFRDWorld](
                new Cpp_EGFRDWorld(
                    deref((<Real3>edge_lengths).thisptr),
                    deref(matrix_sizes.thisptr)))
        elif edge_lengths is None:
            self.thisptr = new shared_ptr[Cpp_EGFRDWorld](new Cpp_EGFRDWorld())
        elif isinstance(edge_lengths, Real3):
            self.thisptr = new shared_ptr[Cpp_EGFRDWorld](
                new Cpp_EGFRDWorld(deref((<Real3>edge_lengths).thisptr)))
        else:
            filename = tostring(edge_lengths)
            self.thisptr = new shared_ptr[Cpp_EGFRDWorld](
                new Cpp_EGFRDWorld(filename))

    def __dealloc__(self):
        # XXX: Here, we release shared pointer,
        #      and if reference count to the EGFRDWorld object,
        #      it will be released automatically.
        del self.thisptr

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

    def set_t(self, Real t):
        """set_t(t)

        Set the value of the time of the world.

        Parameters
        ----------
        t : Real
            The time of the world

        """
        self.thisptr.get().set_t(t)

    def t(self):
        """Return the time of the world."""
        return self.thisptr.get().t()

    def edge_lengths(self):
        """edge_lengths() -> Real3

        Return the edge lengths of the world.

        """
        cdef Cpp_Real3 lengths = self.thisptr.get().edge_lengths()
        return Real3_from_Cpp_Real3(address(lengths))

    def actual_lengths(self):
        """actual_lengths() -> Real3

        Return the actual edge lengths of the world.
        Same as ``edge_lengths``.
        """
        cdef Cpp_Real3 lengths = self.thisptr.get().actual_lengths()
        return Real3_from_Cpp_Real3(address(lengths))

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

    def get_value(self, Species sp):
        """get_value(sp) -> Real

        Return the value (number) corresponding the given Species.

        Parameters
        ----------
        sp : Species
            a species whose value you require

        Returns
        -------
        Real:
            the value

        """
        return self.thisptr.get().get_value(deref(sp.thisptr))

    def get_value_exact(self, Species sp):
        """get_value_exact(sp) -> Real

        Return the value (number) corresponding the given Species.

        Parameters
        ----------
        sp : Species
            a species whose value you require

        Returns
        -------
        Real:
            the value

        """
        return self.thisptr.get().get_value_exact(deref(sp.thisptr))

    def num_particles(self, Species sp = None):
        """num_particles(sp=None) -> Integer

        Return the number of particles.

        Parameters
        ----------
        sp : Species, optional
            The species of particles to count
            If no species is given, return the total number of particles.

        Returns
        -------
            Integer: The number of particles (of the given species)

        """
        if sp is None:
            return self.thisptr.get().num_particles()
        else:
            return self.thisptr.get().num_particles(deref(sp.thisptr))

    def num_particles_exact(self, Species sp):
        """num_particles_exact(sp) -> Integer

        Return the number of particles of a given species.

        Parameters
        ----------
        sp : Species
            The species of particles to count

        Returns
        -------
        Integer:
            The number of particles of a given species

        """
        return self.thisptr.get().num_particles_exact(deref(sp.thisptr))

    def list_particles(self, Species sp = None):
        """list_particles(sp) -> [(ParticleID, Particle)]

        Return the list of particles.

        Parameters
        ----------
        sp : Species, optional
            The species of particles to list up
            If no species is given, return the whole list of particles.

        Returns
        -------
        list:
            The list of particles (of the given species)

        """
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]] particles
        if sp is None:
            particles = self.thisptr.get().list_particles()
        else:
            particles = self.thisptr.get().list_particles(deref(sp.thisptr))

        retval = []
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]].iterator \
            it = particles.begin()
        while it != particles.end():
            retval.append(
                (ParticleID_from_Cpp_ParticleID(
                     <Cpp_ParticleID*>(address(deref(it).first))),
                 Particle_from_Cpp_Particle(
                     <Cpp_Particle*>(address(deref(it).second)))))
            inc(it)
        return retval

    def list_particles_exact(self, Species sp):
        """list_particles_exact(sp) -> [(ParticleID, Particle)]

        Return the list of particles of a given species.

        Parameters
        ----------
        sp : Species
            The species of particles to list up

        Returns
        -------
        list:
            The list of particles of a given species

        """
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]] particles
        particles = self.thisptr.get().list_particles(deref(sp.thisptr))

        retval = []
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]].iterator \
            it = particles.begin()
        while it != particles.end():
            retval.append(
                (ParticleID_from_Cpp_ParticleID(
                     <Cpp_ParticleID*>(address(deref(it).first))),
                 Particle_from_Cpp_Particle(
                     <Cpp_Particle*>(address(deref(it).second)))))
            inc(it)
        return retval

    def has_particle(self, ParticleID pid):
        """has_particle(pid) -> bool

        Check if a particle associated with a given particle id exists.

        Parameters
        ----------
        pid : ParticleID
            A particle id to check

        Returns
        -------
        bool:
            If a particle exists, return True. Otherwise return False

        """
        return self.thisptr.get().has_particle(deref(pid.thisptr))

    def update_particle(self, ParticleID pid, Particle p):
        """update_particle(pid, p) -> bool

        Update a particle.

        Parameters
        ----------
        pid : ParticleID
            A particle id of the particle to update
        p : Particle
            The information to update a particle

        Returns
        -------
        bool:
            True if a new particle was created.

        """
        return self.thisptr.get().update_particle(deref(pid.thisptr), deref(p.thisptr))

    def get_particle(self, ParticleID pid):
        """get_particle(pid) -> (ParticleID, Particle)

        Return the particle associated a given ParticleID.

        Parameters
        ----------
        pid : ParticleID
            An id of the particle you want

        Returns
        -------
        tuple:
            A pair of ParticleID and Particle

        """
        cdef pair[Cpp_ParticleID, Cpp_Particle] \
            pid_particle_pair = self.thisptr.get().get_particle(deref(pid.thisptr))
        return (ParticleID_from_Cpp_ParticleID(address(pid_particle_pair.first)),
                Particle_from_Cpp_Particle(address(pid_particle_pair.second)))

    def remove_particle(self, ParticleID pid):
        """remove_particle(pid)

        Remove the particle associated with a given ParticleID.

        Parameters
        ----------
        pid : ParticleID
            An id of particle to remove

        """
        self.thisptr.get().remove_particle(deref(pid.thisptr))

    def list_particles_within_radius(
        self, Real3 pos, Real radius,
        ParticleID ignore1 = None, ParticleID ignore2 = None):
        """list_particles_within_radius(pos, radius, ignore1=None, ignore2=None) -> [((ParticleID, Particle), Real)]

        Returns a list of pairs of a particle and distance within the given sphere.
        The region is specified with a center position and radius.
        ignore1 and ignore2 will be removed from the list.

        Parameters
        ----------
        pos : Real3
            A center position.
        radius : Real
            A radius.
        ignore1 : ParticleID, optional
            An id ignored.
        ignore2 : ParticleID, optional
            An id ignored.

        Returns
        -------
        list:
            A list of pairs of a particle and its distance from the center position.

        """
        cdef vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]] particles
        if ignore1 is None and ignore2 is None:
            particles = self.thisptr.get().list_particles_within_radius(
                deref(pos.thisptr), radius)
        elif ignore2 is None:
            particles = self.thisptr.get().list_particles_within_radius(
                deref(pos.thisptr), radius, deref(ignore1.thisptr))
        else:
            particles = self.thisptr.get().list_particles_within_radius(
                deref(pos.thisptr), radius,
                deref(ignore1.thisptr), deref(ignore2.thisptr))

        retval = []
        cdef vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]].iterator \
            it = particles.begin()
        while it != particles.end():
            retval.append(
                ((ParticleID_from_Cpp_ParticleID(
                      <Cpp_ParticleID*>(address(deref(it).first.first))),
                  Particle_from_Cpp_Particle(
                      <Cpp_Particle*>(address(deref(it).first.second)))),
                 deref(it).second))
            inc(it)
        return retval

    # def periodic_transpose(self, Real3 pos1, Real3 pos2):
    #     """periodic_transpose(Real3 pos1, Real3 pos2) -> Real3
    #
    #     Return a closest image of pos1 relative to the given position (pos2).
    #
    #     """
    #     cdef Cpp_Real3 newpos = self.thisptr.get().periodic_transpose(
    #         deref(pos1.thisptr), deref(pos2.thisptr))
    #     return Real3_from_Cpp_Real3(address(newpos))

    def apply_boundary(self, Real3 pos):
        """apply_boundary(Real3 pos) -> Real3

        Return a position within the world by applying periodic boundaries
        to the given position.

        """
        cdef Cpp_Real3 newpos = self.thisptr.get().apply_boundary(deref(pos.thisptr))
        return Real3_from_Cpp_Real3(address(newpos))

    # def distance_sq(self, Real3 pos1, Real3 pos2):
    #     """distance_sq(Real3 pos1, Real3 pos2) -> Real
    #
    #     Return a square of the closest distance between the given positions.
    #
    #     """
    #     return self.thisptr.get().distance_sq(deref(pos1.thisptr), deref(pos2.thisptr))

    def distance(self, Real3 pos1, Real3 pos2):
        """distance(Real3 pos1, Real3 pos2) -> Real

        Return the closest distance between the given positions.

        """
        return self.thisptr.get().distance(deref(pos1.thisptr), deref(pos2.thisptr))

    def volume(self):
        """Return the volume of the world."""
        return self.thisptr.get().volume()

    def has_species(self, Species sp):
        """has_species(sp) -> bool

        Check if the given species is in the space or not.

        Parameters
        ----------
        sp : Species
            A species to be found.

        Returns
        -------
        bool:
            True if the species in the space.

        """
        return self.thisptr.get().has_species(deref(sp.thisptr))

    def list_species(self):
        """Return a list of species."""
        cdef vector[Cpp_Species] species = self.thisptr.get().list_species()

        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(
                 Species_from_Cpp_Species(
                     <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    def num_molecules(self, Species sp):
        """num_molecules(sp) -> Integer

        Return the number of molecules.

        Parameters
        ----------
        sp : Species
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

    # def add_species(self, Species sp):
    #     self.thisptr.get().add_species(deref(sp.thisptr))

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

        Remove the molecules.

        Parameters
        ----------
        sp : Species
            a species whose molecules to remove
        num : Integer
            a number of molecules to be removed

        """
        self.thisptr.get().remove_molecules(deref(sp.thisptr), num)

    def save(self, filename):
        """save(filename)

        Save the world to a file.

        Parameters
        ----------
        filename : str
            a filename to save to

        """
        self.thisptr.get().save(tostring(filename))

    def load(self, filename):
        """load(filename)

        Load the world from a file.

        Parameters
        ----------
        filename : str
            a filename to load from

        """
        self.thisptr.get().load(tostring(filename))

    def bind_to(self, m):
        """bind_to(m)

        Bind a model to the world

        Parameters
        ----------
        m : Model
            a model to bind

        """
        self.thisptr.get().bind_to(Cpp_Model_from_Model(m))

    def rng(self):
        """Return a random number generator object."""
        return GSLRandomNumberGenerator_from_Cpp_RandomNumberGenerator(
            self.thisptr.get().rng())

cdef EGFRDWorld EGFRDWorld_from_Cpp_EGFRDWorld(
    shared_ptr[Cpp_EGFRDWorld] w):
    r = EGFRDWorld(Real3(1, 1, 1))
    r.thisptr.swap(w)
    return r

## EGFRDSimulator
#  a python wrapper for Cpp_EGFRDSimulator
cdef class EGFRDSimulator:
    """ A class running the simulation with the egfrd algorithm.

    EGFRDSimulator(m, w)

    """

    def __init__(self, *args):
        """EGFRDSimulator(m, w, bd_dt_factor, dissociation_retry_moves, user_max_shell_size)
        EGFRDSimulator(w, bd_dt_factor, dissociation_retry_moves, user_max_shell_size)

        Constructor.

        Parameters
        ----------
        m : Model
            A model
        w : EGFRDWorld
            A world
        bd_dt_factor : Real
        dissociation_retry_moves : Integer
        user_max_shell_size : Real

        """
        pass

    def __cinit__(self, *args):
        if len(args) == 1:
            self.thisptr = new Cpp_EGFRDSimulator(deref((<EGFRDWorld>args[0]).thisptr))
        elif len(args) == 2:
            if isinstance(args[1], EGFRDWorld):
                self.thisptr = new Cpp_EGFRDSimulator(
                    deref((<EGFRDWorld>args[1]).thisptr),
                    Cpp_Model_from_Model(args[0]))
            else:
                self.thisptr = new Cpp_EGFRDSimulator(
                    deref((<EGFRDWorld>args[0]).thisptr),
                    <Real>args[1])
        elif len(args) == 3:
            if isinstance(args[1], EGFRDWorld):
                self.thisptr = new Cpp_EGFRDSimulator(
                    deref((<EGFRDWorld>args[1]).thisptr),
                    Cpp_Model_from_Model(args[0]),
                    <Real>args[2])
            else:
                self.thisptr = new Cpp_EGFRDSimulator(
                    deref((<EGFRDWorld>args[0]).thisptr),
                    <Real>args[1], <Integer>args[2])
        elif len(args) == 4:
            if isinstance(args[1], EGFRDWorld):
                self.thisptr = new Cpp_EGFRDSimulator(
                    deref((<EGFRDWorld>args[1]).thisptr),
                    Cpp_Model_from_Model(args[0]),
                    <Real>args[2], <Integer>args[3])
            else:
                self.thisptr = new Cpp_EGFRDSimulator(
                    deref((<EGFRDWorld>args[0]).thisptr),
                    <Real>args[1], <Integer>args[2], <Real>args[3])
        elif len(args) == 5:
            self.thisptr = new Cpp_EGFRDSimulator(
                deref((<EGFRDWorld>args[1]).thisptr),
                Cpp_Model_from_Model(args[0]),
                <Real>args[2], <Integer>args[3], <Real>args[4])
        else:
            raise ValueError(
                "The invalid number of arguments was given [{}].".format(len(args)))

    def __dealloc__(self):
        del self.thisptr

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.num_steps()

    def step(self, upto = None):
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
            return self.thisptr.step(<Real> upto)

    def t(self):
        """Return the time."""
        return self.thisptr.t()

    def dt(self):
        """Return the step interval."""
        return self.thisptr.dt()

    def next_time(self):
        """Return the scheduled time for the next step."""
        return self.thisptr.next_time()

    def check_reaction(self):
        """Return if any reaction occurred at the last step, or not."""
        return self.thisptr.check_reaction()

    def last_reactions(self):
        """last_reactions() -> [(ReactionRule, ReactionInfo)]

        Return reactions occuring at the last step.

        Returns
        -------
        list:
            the list of reaction rules and infos.

        """
        cdef vector[pair[Cpp_ReactionRule, Cpp_ReactionInfo]] reactions = self.thisptr.last_reactions()
        cdef vector[pair[Cpp_ReactionRule, Cpp_ReactionInfo]].iterator it = reactions.begin()
        retval = []
        while it != reactions.end():
            retval.append((
                ReactionRule_from_Cpp_ReactionRule(
                    <Cpp_ReactionRule*>(address(deref(it).first))),
                ReactionInfo_from_Cpp_ReactionInfo(
                    <Cpp_ReactionInfo*>(address(deref(it).second)))))
            inc(it)
        return retval


    def set_t(self, Real new_t):
        """set_t(t)

        Set the current time.

        Parameters
        ----------
        t : Real
            A current time.

        """
        self.thisptr.set_t(new_t)

    def set_dt(self, Real dt):
        """set_dt(dt)

        Set a step interval.

        Parameters
        ----------
        dt : Real
            A step interval

        """
        self.thisptr.set_dt(dt)

    def initialize(self):
        """Initialize the simulator."""
        self.thisptr.initialize()

    def model(self):
        """Return the model bound."""
        return Model_from_Cpp_Model(self.thisptr.model())

    def world(self):
        """Return the world bound."""
        return EGFRDWorld_from_Cpp_EGFRDWorld(self.thisptr.world())

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

    def set_paranoiac(self, val):
        self.thisptr.set_paranoiac(<bool>val)

cdef EGFRDSimulator EGFRDSimulator_from_Cpp_EGFRDSimulator(Cpp_EGFRDSimulator* s):
    r = EGFRDSimulator(
        Model_from_Cpp_Model(s.model()), EGFRDWorld_from_Cpp_EGFRDWorld(s.world()))
    del r.thisptr
    r.thisptr = s
    return r

## EGFRDFactory
#  a python wrapper for Cpp_EGFRDFactory
cdef class EGFRDFactory:
    """ A factory class creating a EGFRDWorld instance and a EGFRDSimulator instance.

    EGFRDFactory(Integer3 matrix_sizes=None, Real bd_dt_factor=None,
                 Integer dissociation_retry_moves=None, Real user_max_shell_size=None)

    """

    def __init__(self, Integer3 matrix_sizes=None, bd_dt_factor=None,
                 dissociation_retry_moves=None, user_max_shell_size=None):
        """Constructor.

        Parameters
        ----------
        matrix_sizes : Integer3, optional
            A size of a cell matrix.
            The number of cells must be larger than 3, in principle.
        bd_dt_factor : Real, optioanl
            A rescaling factor for the step interval
            of BD propagation in a Multi domain.
        dissociation_retry_moves : Integer, optional
            A number of trials for placing a new product when it's failed
            because of the overlap.
        user_max_shell_size : Real, optional
            A custom max shell size.

        """
        pass

    def __cinit__(self, Integer3 matrix_sizes=None, bd_dt_factor=None,
                  dissociation_retry_moves=None, user_max_shell_size=None):
        self.thisptr = new Cpp_EGFRDFactory(
            Cpp_EGFRDFactory.default_matrix_sizes() if matrix_sizes is None else deref(matrix_sizes.thisptr),
            Cpp_EGFRDFactory.default_bd_dt_factor() if bd_dt_factor is None else <Real>bd_dt_factor,
            Cpp_EGFRDFactory.default_dissociation_retry_moves() if dissociation_retry_moves is None else <Integer>dissociation_retry_moves,
            Cpp_EGFRDFactory.default_user_max_shell_size() if user_max_shell_size is None else <Real>user_max_shell_size)

    def __dealloc__(self):
        del self.thisptr

    def rng(self, GSLRandomNumberGenerator rng):
        """rng(GSLRandomNumberGenerator) -> EGFRDFactory

        Set a random number generator, and return self.

        """
        cdef Cpp_EGFRDFactory *ptr = self.thisptr.rng_ptr(deref(rng.thisptr))
        assert ptr == self.thisptr
        return self

    def create_world(self, arg1=None):
        """create_world(arg1=None) -> EGFRDWorld

        Return a EGFRDWorld instance.

        Parameters
        ----------
        arg1 : Real3
            The lengths of edges of a EGFRDWorld created

        or

        arg1 : str
            The path of a HDF5 file for EGFRDWorld

        Returns
        -------
        EGFRDWorld:
            The created world

        """
        if arg1 is None:
            return EGFRDWorld_from_Cpp_EGFRDWorld(
                shared_ptr[Cpp_EGFRDWorld](self.thisptr.create_world()))
        elif isinstance(arg1, Real3):
            return EGFRDWorld_from_Cpp_EGFRDWorld(
                shared_ptr[Cpp_EGFRDWorld](
                    self.thisptr.create_world(deref((<Real3>arg1).thisptr))))
        elif isinstance(arg1, str):
            return EGFRDWorld_from_Cpp_EGFRDWorld(
                shared_ptr[Cpp_EGFRDWorld](self.thisptr.create_world(<string>(arg1))))
        else:
            return EGFRDWorld_from_Cpp_EGFRDWorld(
                shared_ptr[Cpp_EGFRDWorld](self.thisptr.create_world(
                    Cpp_Model_from_Model(arg1))))

    def create_simulator(self, arg1, EGFRDWorld arg2=None):
        """create_simulator(arg1, arg2) -> EGFRDSimulator

        Return a EGFRDSimulator instance.

        Parameters
        ----------
        arg1 : EGFRDWorld
            A world

        or

        arg1 : Model
            A simulation model
        arg2 : EGFRDWorld
            A world

        Returns
        -------
        EGFRDSimulator:
            The created simulator

        """
        if arg2 is None:
            return EGFRDSimulator_from_Cpp_EGFRDSimulator(
                self.thisptr.create_simulator(deref((<EGFRDWorld>arg1).thisptr)))
        else:
            return EGFRDSimulator_from_Cpp_EGFRDSimulator(
                self.thisptr.create_simulator(
                    Cpp_Model_from_Model(arg1), deref(arg2.thisptr)))

## BDSimulator
#  a python wrapper for Cpp_BDSimulator
cdef class BDSimulator:
    """ A class running the simulation with the bd algorithm.

    BDSimulator(m, w)

    """

    def __init__(self, *args):
        """BDSimulator(m, w, bd_dt_factor, dissociation_retry_moves)
        BDSimulator(w, bd_dt_factor, dissociation_retry_moves)

        Constructor.

        Parameters
        ----------
        m : Model
            A model
        w : EGFRDWorld
            A world
        bd_dt_factor : Real
        dissociation_retry_moves : Integer

        """
        pass

    def __cinit__(self, *args):
        if len(args) == 1:
            self.thisptr = new Cpp_BDSimulator(deref((<EGFRDWorld>args[0]).thisptr))
        elif len(args) == 2:
            if isinstance(args[1], EGFRDWorld):
                self.thisptr = new Cpp_BDSimulator(
                    deref((<EGFRDWorld>args[1]).thisptr),
                    Cpp_Model_from_Model(args[0]))
            else:
                self.thisptr = new Cpp_BDSimulator(
                    deref((<EGFRDWorld>args[0]).thisptr),
                    <Real>args[1])
        elif len(args) == 3:
            if isinstance(args[1], EGFRDWorld):
                self.thisptr = new Cpp_BDSimulator(
                    deref((<EGFRDWorld>args[1]).thisptr),
                    Cpp_Model_from_Model(args[0]),
                    <Real>args[2])
            else:
                self.thisptr = new Cpp_BDSimulator(
                    deref((<EGFRDWorld>args[0]).thisptr),
                    <Real>args[1], <Integer>args[2])
        elif len(args) == 4:
            self.thisptr = new Cpp_BDSimulator(
                deref((<EGFRDWorld>args[1]).thisptr),
                Cpp_Model_from_Model(args[0]),
                <Real>args[2], <Integer>args[3])
        else:
            raise ValueError(
                "The invalid number of arguments was given [{}].".format(len(args)))

    def __dealloc__(self):
        del self.thisptr

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.num_steps()

    def step(self, upto = None):
        """step(upto=None) -> bool

        Step the simulation.

        Parameters
        ----------
        upto : Real, optional
            The time which to step the simulation up to

        Returns
        -------
        bool:
            True if the simulation did not reach the given time.
            When upto is not given, nothing will be returned.

        """
        if upto is None:
            self.thisptr.step()
        else:
            return self.thisptr.step(<Real> upto)

    def t(self):
        """Return the time."""
        return self.thisptr.t()

    def dt(self):
        """Return the step interval."""
        return self.thisptr.dt()

    def next_time(self):
        """Return the scheduled time for the next step."""
        return self.thisptr.next_time()

    def check_reaction(self):
        """Return if any reaction occurred at the last step, or not."""
        return self.thisptr.check_reaction()

    def last_reactions(self):
        """last_reactions() -> [(ReactionRule, ReactionInfo)]

        Return reactions occuring at the last step.

        Returns
        -------
        list:
            The list of reaction rules and infos.

        """
        cdef vector[pair[Cpp_ReactionRule, Cpp_ReactionInfo]] reactions = self.thisptr.last_reactions()
        cdef vector[pair[Cpp_ReactionRule, Cpp_ReactionInfo]].iterator it = reactions.begin()
        retval = []
        while it != reactions.end():
            retval.append((
                ReactionRule_from_Cpp_ReactionRule(
                    <Cpp_ReactionRule*>(address(deref(it).first))),
                ReactionInfo_from_Cpp_ReactionInfo(
                    <Cpp_ReactionInfo*>(address(deref(it).second)))))
            inc(it)
        return retval

    def set_t(self, Real new_t):
        """set_t(t)

        Set the current time.

        Parameters
        ----------
        t : Real
            a current time.

        """
        self.thisptr.set_t(new_t)

    def set_dt(self, Real dt):
        """set_dt(dt)

        Set a step interval.

        Parameters
        ----------
        dt : Real
            a step interval

        """
        self.thisptr.set_dt(dt)

    def initialize(self):
        """Initialize the simulator."""
        self.thisptr.initialize()

    def model(self):
        """Return the model bound."""
        return Model_from_Cpp_Model(self.thisptr.model())

    def world(self):
        """Return the world bound."""
        return EGFRDWorld_from_Cpp_EGFRDWorld(self.thisptr.world())

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

    def dt_factor(self):
        return self.thisptr.dt_factor()

    def add_potential(self, Species sp, shape, threshold=None):
        """add_potential(sp, arg)

        Set the potential for each Species.

        Parameters
        ----------
        sp : Species
            a species of molecules to add
        shape : Real or Shape
            a radius or a shape for the potential
        threshold : Real, optional
            a threshold
            default is None.

        """
        if threshold is None:
            if isinstance(shape, numbers.Number):
                self.thisptr.add_potential(deref(sp.thisptr), <Real>shape)
            else:
                self.thisptr.add_potential(
                    deref(sp.thisptr), deref((<Shape>(shape.as_base())).thisptr))
        else:
            self.thisptr.add_potential(
                deref(sp.thisptr), deref((<Shape>(shape.as_base())).thisptr), <Real>threshold)


cdef BDSimulator BDSimulator_from_Cpp_BDSimulator(Cpp_BDSimulator* s):
    r = BDSimulator(
        Model_from_Cpp_Model(s.model()), EGFRDWorld_from_Cpp_EGFRDWorld(s.world()))
    del r.thisptr
    r.thisptr = s
    return r

## BDFactory
#  a python wrapper for Cpp_BDFactory
cdef class BDFactory:
    """ A factory class creating a BDWorld instance and a BDSimulator instance.

    BDFactory(Integer3 matrix_sizes=None, Real bd_dt_factor=None,
              Integer dissociation_retry_moves=None)

    """

    def __init__(self, Integer3 matrix_sizes=None, bd_dt_factor=None,
                 dissociation_retry_moves=None):
        """Constructor.

        Parameters
        ----------
        matrix_sizes : Integer3, optional
            A size of a cell matrix.
            The number of cells must be larger than 3, in principle.
        bd_dt_factor : Real, optioanl
            A rescaling factor for the step interval
            of BD propagation in a Multi domain.
        dissociation_retry_moves : Integer, optional
            A number of trials for placing a new product when it's failed
            because of the overlap.

        """
        pass

    def __cinit__(self, Integer3 matrix_sizes=None, bd_dt_factor=None,
                  dissociation_retry_moves=None):
        self.thisptr = new Cpp_BDFactory(
            Cpp_BDFactory.default_matrix_sizes() if matrix_sizes is None else deref(matrix_sizes.thisptr),
            Cpp_BDFactory.default_bd_dt_factor() if bd_dt_factor is None else <Real>bd_dt_factor,
            Cpp_BDFactory.default_dissociation_retry_moves() if dissociation_retry_moves is None else <Integer>dissociation_retry_moves)

    def __dealloc__(self):
        del self.thisptr

    def rng(self, GSLRandomNumberGenerator rng):
        """rng(GSLRandomNumberGenerator) -> BDFactory

        Set a random number generator, and return self.

        """
        cdef Cpp_BDFactory *ptr = self.thisptr.rng_ptr(deref(rng.thisptr))
        assert ptr == self.thisptr
        return self

    def create_world(self, arg1):
        """create_world(arg1=None) -> EGFRDWorld

        Return a EGFRDWorld instance.

        Parameters
        ----------
        arg1 : Real3
            The lengths of edges of a EGFRDWorld created

        or

        arg1 : str
            The path of a HDF5 file for EGFRDWorld

        Returns
        -------
        EGFRDWorld:
            The created world

        """
        if isinstance(arg1, Real3):
            return EGFRDWorld_from_Cpp_EGFRDWorld(
                shared_ptr[Cpp_EGFRDWorld](
                    self.thisptr.create_world(deref((<Real3>arg1).thisptr))))
        elif isinstance(arg1, str):
            return EGFRDWorld_from_Cpp_EGFRDWorld(
                shared_ptr[Cpp_EGFRDWorld](self.thisptr.create_world(<string>(arg1))))
        else:
            return EGFRDWorld_from_Cpp_EGFRDWorld(
                shared_ptr[Cpp_EGFRDWorld](self.thisptr.create_world(
                    Cpp_Model_from_Model(arg1))))

    def create_simulator(self, arg1, EGFRDWorld arg2=None):
        """create_simulator(arg1, arg2) -> BDSimulator

        Return a BDSimulator instance.

        Parameters
        ----------
        arg1 : BDWorld
            A world

        or

        arg1 : Model
            A simulation model
        arg2 : BDWorld
            A world

        Returns
        -------
        BDSimulator:
            The created simulator

        """
        if arg2 is None:
            return BDSimulator_from_Cpp_BDSimulator(
                self.thisptr.create_simulator(deref((<EGFRDWorld>arg1).thisptr)))
        else:
            return BDSimulator_from_Cpp_BDSimulator(
                self.thisptr.create_simulator(
                    Cpp_Model_from_Model(arg1), deref(arg2.thisptr)))

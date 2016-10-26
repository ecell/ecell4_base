import collections
from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address
from libcpp.string cimport string
from libcpp.vector cimport vector

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *


## ReactionInfo
cdef class ReactionInfo:
    """A class stores detailed information about a reaction in bd.

    ReactionInfo(t, reactants, products)

    """

    def __init__(self, Real t, reactants, products):
        """Constructor.

        Args:
          t (Real): A time when a reaction occurred
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

## BDWorld
#  a python wrapper for Cpp_BDWorld
cdef class BDWorld:
    """A class containing the properties of the bd world.

    BDWorld(edge_lengths=None, matrix_sizes=None, GSLRandomNumberGenerator rng=None)

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

        if edge_lengths is None:
            self.thisptr = new shared_ptr[Cpp_BDWorld](new Cpp_BDWorld())
        elif matrix_sizes is None:
            if isinstance(edge_lengths, Real3):
                self.thisptr = new shared_ptr[Cpp_BDWorld](
                    new Cpp_BDWorld(deref((<Real3>edge_lengths).thisptr)))
            else:
                filename = tostring(edge_lengths)
                self.thisptr = new shared_ptr[Cpp_BDWorld](new Cpp_BDWorld(filename))
        elif rng is None:
            self.thisptr = new shared_ptr[Cpp_BDWorld](
                new Cpp_BDWorld(deref((<Real3>edge_lengths).thisptr),
                    deref(matrix_sizes.thisptr)))
        else:
            self.thisptr = new shared_ptr[Cpp_BDWorld](
                new Cpp_BDWorld(deref((<Real3>edge_lengths).thisptr),
                    deref(matrix_sizes.thisptr), deref(rng.thisptr)))

    def __dealloc__(self):
        # XXX: Here, we release shared pointer,
        #      and if reference count to the BDWorld object,
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
        Integer:
            The number of particles (of the given species)

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

    def periodic_transpose(self, Real3 pos1, Real3 pos2):
        """periodic_transpose(Real3 pos1, Real3 pos2) -> Real3

        Return a closest image of pos1 relative to the given position (pos2).

        """
        cdef Cpp_Real3 newpos = self.thisptr.get().periodic_transpose(
            deref(pos1.thisptr), deref(pos2.thisptr))
        return Real3_from_Cpp_Real3(address(newpos))

    def apply_boundary(self, Real3 pos):
        """apply_boundary(Real3 pos) -> Real3

        Return a position within the world by applying periodic boundaries
        to the given position.

        """
        cdef Cpp_Real3 newpos = self.thisptr.get().apply_boundary(deref(pos.thisptr))
        return Real3_from_Cpp_Real3(address(newpos))

    def distance_sq(self, Real3 pos1, Real3 pos2):
        """distance_sq(Real3 pos1, Real3 pos2) -> Real

        Return a square of the closest distance between the given positions.

        """
        return self.thisptr.get().distance_sq(deref(pos1.thisptr), deref(pos2.thisptr))

    def distance(self, Real3 pos1, Real3 pos2):
        """distance(Real3 pos1, Real3 pos2) -> Real

        Return the closest distance between the given positions.

        """
        return self.thisptr.get().distance(deref(pos1.thisptr), deref(pos2.thisptr))

    def volume(self):
        """Return the volume of the world."""
        return self.thisptr.get().volume()

    # def has_species(self, Species sp):
    #     """has_species(sp) -> bool
    #
    #     Check if the given species is in the space or not.
    #
    #     Args:
    #         sp (Species): A species to be found.
    #
    #     Returns:
    #         bool: True if the species in the space.
    #
    #     """
    #     return self.thisptr.get().has_species(deref(sp.thisptr))

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

    def num_molecules(self, Species sp):
        """num_molecules(sp) -> Integer

        Return the number of molecules.

        Parameters
        ----------
        sp : Species
            A species whose molecules you count

        Returns
        -------
        Integer:
            The number of molecules (of a given species)

        """
        # if sp is None:
        #     return self.thisptr.get().num_molecules()
        # else:
        #     return self.thisptr.get().num_molecules(deref(sp.thisptr))
        return self.thisptr.get().num_molecules(deref(sp.thisptr))

    def num_molecules_exact(self, Species sp):
        """num_molecules_exact(sp) -> Integer

        Return the number of molecules of a given species.

        Parameters
        ----------
        sp : Species
            A species whose molecules you count

        Returns
        -------
        Integer:
            The number of molecules of a given species

        """
        return self.thisptr.get().num_molecules_exact(deref(sp.thisptr))

    # def add_species(self, Species sp):
    #     self.thisptr.get().add_species(deref(sp.thisptr))

    # def add_molecules(self, Species sp, Integer num):
    #     self.thisptr.get().add_molecules(deref(sp.thisptr), num)

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

    def as_base(self):
        """Return self as a base class. Only for developmental use."""
        retval = Space()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Space](
            <shared_ptr[Cpp_Space]>deref(self.thisptr))
        return retval

cdef BDWorld BDWorld_from_Cpp_BDWorld(
    shared_ptr[Cpp_BDWorld] w):
    r = BDWorld(Real3(1, 1, 1))
    r.thisptr.swap(w)
    return r

## BDSimulator
#  a python wrapper for Cpp_BDSimulator
cdef class BDSimulator:
    """ A class running the simulation with the bd algorithm.

    BDSimulator(m, w, bd_dt_factor)

    """

    def __init__(self, m, BDWorld w=None, bd_dt_factor=None):
        """BDSimulator(m, w, bd_dt_factor)
        BDSimulator(w, bd_dt_factor)

        Constructor.

        Parameters
        ----------
        m : Model
            A model
        w : BDWorld
            A world
        bd_dt_factor : Real

        """
        pass

    def __cinit__(self, m, BDWorld w=None, bd_dt_factor=None):
        if w is None:
            if bd_dt_factor is None:
                self.thisptr = new Cpp_BDSimulator(
                    deref((<BDWorld>m).thisptr))
            else:
                self.thisptr = new Cpp_BDSimulator(
                    deref((<BDWorld>m).thisptr), <Real>bd_dt_factor)
        else:
            if bd_dt_factor is None:
                self.thisptr = new Cpp_BDSimulator(
                    Cpp_Model_from_Model(m), deref(w.thisptr))
            else:
                self.thisptr = new Cpp_BDSimulator(
                    Cpp_Model_from_Model(m), deref(w.thisptr),
                    <Real>bd_dt_factor)

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
            return self.thisptr.step(upto)

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

    def set_dt(self, Real& dt):
        """set_dt(dt)

        Set a step interval.

        Parameters
        ----------
        dt : Real
            a step interval

        """
        self.thisptr.set_dt(dt)

    def next_time(self):
        """Return the scheduled time for the next step."""
        return self.thisptr.next_time()

    def initialize(self):
        """Initialize the simulator."""
        self.thisptr.initialize()

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

    def model(self):
        """Return the model bound."""
        return Model_from_Cpp_Model(self.thisptr.model())

    def world(self):
        """Return the world bound."""
        return BDWorld_from_Cpp_BDWorld(self.thisptr.world())

    def run(self, Real duration, observers=None):
        """run(duration, observers)

        Run the simulation.

        Parameters
        ----------
        duration : Real
            a duration for running a simulation.
            A simulation is expected to be stopped at ``t() + duration``.
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

cdef BDSimulator BDSimulator_from_Cpp_BDSimulator(Cpp_BDSimulator* s):
    r = BDSimulator(
        Model_from_Cpp_Model(s.model()), BDWorld_from_Cpp_BDWorld(s.world()))
    del r.thisptr
    r.thisptr = s
    return r

## BDFactory
#  a python wrapper for Cpp_BDFactory
cdef class BDFactory:
    """ A factory class creating a BDWorld instance and a BDSimulator instance.

    BDFactory(Integer3 matrix_sizes=None, Real bd_dt_factor=None)

    """

    def __init__(self, Integer3 matrix_sizes=None, bd_dt_factor=None):
        """Constructor.

        Parameters
        ----------
        matrix_sizes : Integer3, optional
            A size of a cell matrix.
            The number of cells must be larger than 3, in principle.
        bd_dt_factor : Real

        """
        pass

    def __cinit__(self, Integer3 matrix_sizes=None, bd_dt_factor=None):
        self.thisptr = new Cpp_BDFactory(
            Cpp_BDFactory.default_matrix_sizes() if matrix_sizes is None else deref(matrix_sizes.thisptr),
            Cpp_BDFactory.default_bd_dt_factor() if bd_dt_factor is None else <Real>bd_dt_factor)

    def __dealloc__(self):
        del self.thisptr

    def rng(self, GSLRandomNumberGenerator rng):
        """rng(GSLRandomNumberGenerator) -> BDFactory

        Set a random number generator, and return self.

        """
        cdef Cpp_BDFactory *ptr = self.thisptr.rng_ptr(deref(rng.thisptr))
        assert ptr == self.thisptr
        return self

    def create_world(self, arg1=None):
        """create_world(arg1=None) -> BDWorld

        Return a ``BDWorld`` instance.

        Parameters
        ----------
        arg1 : Real3
            The lengths of edges of a ``BDWorld`` created

        or

        arg1 : str
            The path of a HDF5 file for ``BDWorld``

        Returns
        -------
        BDWorld:
            The created world

        """
        if arg1 is None:
            return BDWorld_from_Cpp_BDWorld(
                shared_ptr[Cpp_BDWorld](self.thisptr.create_world()))
        elif isinstance(arg1, Real3):
            return BDWorld_from_Cpp_BDWorld(
                shared_ptr[Cpp_BDWorld](
                    self.thisptr.create_world(deref((<Real3>arg1).thisptr))))
        elif isinstance(arg1, str):
            return BDWorld_from_Cpp_BDWorld(
                shared_ptr[Cpp_BDWorld](self.thisptr.create_world(<string>(arg1))))
        else:
            return BDWorld_from_Cpp_BDWorld(
                shared_ptr[Cpp_BDWorld](self.thisptr.create_world(
                    Cpp_Model_from_Model(arg1))))

    def create_simulator(self, arg1, BDWorld arg2=None):
        """create_simulator(arg1, arg2=None) -> BDSimulator

        Return a ``BDSimulator`` instance.

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
                self.thisptr.create_simulator(deref((<BDWorld>arg1).thisptr)))
        else:
            return BDSimulator_from_Cpp_BDSimulator(
                self.thisptr.create_simulator(
                    Cpp_Model_from_Model(arg1), deref(arg2.thisptr)))

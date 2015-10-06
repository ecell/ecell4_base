import collections
from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address
from libcpp.string cimport string
from libcpp.vector cimport vector

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *


## BDWorld
#  a python wrapper for Cpp_BDWorld
cdef class BDWorld:
    """A class containing the properties of the bd world.

    BDWorld(edge_lengths=None, matrix_sizes=None, GSLRandomNumberGenerator rng=None)

    """

    def __init__(self, edge_lengths = None, Integer3 matrix_sizes = None,
                 GSLRandomNumberGenerator rng = None):
        """Constructor.

        Args:
            edge_lengths (Real3, optional): A size of the World.
            matrix_sizes (Integer3, optional): A size of a cell matrix.
                The number of cells must be larger than 3, in principle.
            rng (GSLRandomNumberGenerator, optional): A random number generator.

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

        Args:
            arg1 (Particle): A particle to be placed.

            or

            arg1 (Species): A species of a particle
            arg2 (Real3): A coordinate to place a particle

        Returns:
            tuple: A pair of ParticleID and Particle of a new particle

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

        Args:
            t (Real): The time of the world

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

    def num_particles(self, Species sp = None):
        """num_particles(sp=None) -> Integer

        Return the number of particles.

        Args:
            sp (Species, optional): The species of particles to count
                If no species is given, return the total number of particles.

        Returns:
            Integer: The number of particles (of the given species)

        """
        if sp is None:
            return self.thisptr.get().num_particles()
        else:
            return self.thisptr.get().num_particles(deref(sp.thisptr))

    def num_particles_exact(self, Species sp):
        """num_particles_exact(sp) -> Integer

        Return the number of particles of a given species.

        Args:
            sp (Species): The species of particles to count

        Returns:
            Integer: The number of particles of a given species

        """
        return self.thisptr.get().num_particles_exact(deref(sp.thisptr))

    def list_particles(self, Species sp = None):
        """list_particles(sp) -> [(ParticleID, Particle)]

        Return the list of particles.

        Args:
            sp (Species, optional): The species of particles to list up
                If no species is given, return the whole list of particles.

        Returns:
            list: The list of particles (of the given species)

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

        Args:
            sp (Species): The species of particles to list up

        Returns:
            list: The list of particles of a given species

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

        Args:
            pid (ParticleID): A particle id to check

        Returns:
            bool: If a particle exists, return True. Otherwise return False

        """
        return self.thisptr.get().has_particle(deref(pid.thisptr))

    def update_particle(self, ParticleID pid, Particle p):
        """update_particle(pid, p) -> bool

        Update a particle.

        Args:
            pid (ParticleID): A particle id of the particle to update
            p (Particle): The information to update a particle

        Returns:
            bool: True if a new particle was created.

        """
        return self.thisptr.get().update_particle(deref(pid.thisptr), deref(p.thisptr))

    def get_particle(self, ParticleID pid):
        """get_particle(pid) -> (ParticleID, Particle)

        Return the particle associated a given ParticleID.

        Args:
            pid (ParticleID): An id of the particle you want

        Returns:
            tuple: A pair of ParticleID and Particle

        """
        cdef pair[Cpp_ParticleID, Cpp_Particle] \
            pid_particle_pair = self.thisptr.get().get_particle(deref(pid.thisptr))
        return (ParticleID_from_Cpp_ParticleID(address(pid_particle_pair.first)),
                Particle_from_Cpp_Particle(address(pid_particle_pair.second)))

    def remove_particle(self, ParticleID pid):
        """remove_particle(pid)

        Remove the particle associated with a given ParticleID.

        Args:
            pid (ParticleID): An id of particle to remove

        """
        self.thisptr.get().remove_particle(deref(pid.thisptr))

    def list_particles_within_radius(
        self, Real3 pos, Real radius,
        ParticleID ignore1 = None, ParticleID ignore2 = None):
        """list_particles_within_radius(pos, radius, ignore1=None, ignore2=None)
            -> [((ParticleID, Particle), Real)]

        Returns a list of pairs of a particle and distance within the given sphere.
        The region is specified with a center position and radius.
        ignore1 and ignore2 will be removed from the list.

        Args:
            pos (Real3): A center position.
            radius (Real): A radius.
            ignore1 (ParticleID, optional): An id ignored.
            ignore2 (ParticleID, optional): An id ignored.

        Returns:
            list: A list of pairs of a particle and its distance
                from the center position.

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
    #     return self.thisptr.get().has_species(deref(sp.thisptr))

    def num_molecules(self, Species sp):
        """num_molecules(sp) -> Integer

        Return the number of molecules.

        Args:
            sp (Species, optional): a species whose molecules you count

        Returns:
            Integer: the number of molecules (of a given species)

        """
        # if sp is None:
        #     return self.thisptr.get().num_molecules()
        # else:
        #     return self.thisptr.get().num_molecules(deref(sp.thisptr))
        return self.thisptr.get().num_molecules(deref(sp.thisptr))

    def num_molecules_exact(self, Species sp):
        """num_molecules_exact(sp) -> Integer

        Return the number of molecules of a given species.

        Args:
            sp (Species): a species whose molecules you count

        Returns:
            Integer: the number of molecules of a given species

        """
        return self.thisptr.get().num_molecules_exact(deref(sp.thisptr))

    # def add_species(self, Species sp):
    #     self.thisptr.get().add_species(deref(sp.thisptr))

    # def add_molecules(self, Species sp, Integer num):
    #     self.thisptr.get().add_molecules(deref(sp.thisptr), num)

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

        Remove the molecules.

        Args:
            sp (Species): a species whose molecules to remove
            num (Integer): a number of molecules to be removed

        """
        self.thisptr.get().remove_molecules(deref(sp.thisptr), num)

    def save(self, filename):
        """save(filename)

        Save the world to a file.

        Args:
            filename (str): a filename to save to

        """
        self.thisptr.get().save(tostring(filename))

    def load(self, filename):
        """load(filename)

        Load the world from a file.

        Args:
            filename (str): a filename to load from

        """
        self.thisptr.get().load(tostring(filename))

    def bind_to(self, m):
        """bind_to(m)

        Bind a model to the world

        Args:
            m (Model): a model to bind

        """
        self.thisptr.get().bind_to(deref(Cpp_Model_from_Model(m)))

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
    """ A class running the simulation with the gillespie algorithm.

    GillespieSimulator(m, w)

    """


    def __init__(self, m, BDWorld w=None):
        """BDSimulator(m, w)
        BDSimulator(w)

        Constructor.

        Args:
            m (Model): A model
            w (BDWorld): A world

        """
        pass

    def __cinit__(self, m, BDWorld w=None):
        if w is None:
            self.thisptr = new Cpp_BDSimulator(
                deref((<BDWorld>m).thisptr))
        else:
            self.thisptr = new Cpp_BDSimulator(
                deref(Cpp_Model_from_Model(m)), deref(w.thisptr))

    def __dealloc__(self):
        del self.thisptr

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.num_steps()

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
            self.thisptr.step()
        else:
            return self.thisptr.step(upto)

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

    def set_dt(self, Real& dt):
        """set_dt(dt)

        Set a step interval.

        Args:
            dt (Real): a step interval

        """
        self.thisptr.set_dt(dt)

    def next_time(self):
        """Return the scheduled time for the next step."""
        return self.thisptr.next_time()

    def initialize(self):
        """Initialize the simulator."""
        self.thisptr.initialize()

    def last_reactions(self):
        """last_reactions() -> [(ReactionRule, ReactionInfo)]

        Return reactions occuring at the last step.

        Returns:
            list: the list of reaction rules and infos.

        """
        cdef vector[Cpp_ReactionRule] reactions = self.thisptr.last_reactions()
        cdef vector[Cpp_ReactionRule].iterator it = reactions.begin()
        retval = []
        while it != reactions.end():
            retval.append(ReactionRule_from_Cpp_ReactionRule(
                <Cpp_ReactionRule*>(address(deref(it)))))
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

    BDFactory(Integer3 matrix_sizes=None, GSLRandomNumberGenerator rng=None)

    """

    def __init__(self, Integer3 matrix_sizes=None, GSLRandomNumberGenerator rng=None):
        """Constructor.

        Args:
            matrix_sizes (Integer3, optional): A size of a cell matrix.
                The number of cells must be larger than 3, in principle.
            rng (GSLRandomNumberGenerator, optional): a random number generator.

        """
        pass

    def __cinit__(self, Integer3 matrix_sizes=None, GSLRandomNumberGenerator rng=None):
        if matrix_sizes is None:
            self.thisptr = new Cpp_BDFactory()
        elif rng is None:
            self.thisptr = new Cpp_BDFactory(deref(matrix_sizes.thisptr))
        else:
            self.thisptr = new Cpp_BDFactory(deref(matrix_sizes.thisptr), deref(rng.thisptr))

    def __dealloc__(self):
        del self.thisptr

    def create_world(self, arg1):
        """create_world(arg1=None) -> BDWorld

        Return a BDWorld instance.

        Args:
            arg1 (Real3): The lengths of edges of a BDWorld created

            or

            arg1 (str): The path of a HDF5 file for BDWorld

        Returns:
            BDWorld: the created world

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
                    deref(Cpp_Model_from_Model(arg1)))))

    def create_simulator(self, arg1, BDWorld arg2=None):
        """create_simulator(arg1, arg2) -> BDSimulator

        Return a BDSimulator instance.

        Args:
            arg1 (BDWorld): a world

            or

            arg1 (Model): a simulation model
            arg2 (BDWorld): a world

        Returns:
            BDSimulator: the created simulator

        """
        if arg2 is None:
            return BDSimulator_from_Cpp_BDSimulator(
                self.thisptr.create_simulator(deref((<BDWorld>arg1).thisptr)))
        else:
            return BDSimulator_from_Cpp_BDSimulator(
                self.thisptr.create_simulator(
                    deref(Cpp_Model_from_Model(arg1)), deref(arg2.thisptr)))

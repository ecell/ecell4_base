import collections
from cython cimport address
from cython.operator cimport dereference as deref, preincrement as inc
from ecell4.core cimport *


## GillespieWorld
#  a python wrapper for Cpp_GillespieWorld
cdef class GillespieWorld:
    """A class containing the properties of the gillespie world.

    GillespieWorld(edge_lengths=None, GSLRandomNumberGenerator rng=None)

    """

    def __init__(self, edge_lengths = None, GSLRandomNumberGenerator rng = None):
        """Constructor.

        Args:
            edge_lengths (Real3, optional): A size of the World.
            rng (GSLRandomNumberGenerator, optional): A random number generator.

        """
        pass

    def __cinit__(self, edge_lengths = None, GSLRandomNumberGenerator rng = None):
        cdef string filename

        if edge_lengths is None:
            self.thisptr = new shared_ptr[Cpp_GillespieWorld](new Cpp_GillespieWorld())
        elif rng is None:
            if isinstance(edge_lengths, Real3):
                self.thisptr = new shared_ptr[Cpp_GillespieWorld](
                    new Cpp_GillespieWorld(deref((<Real3>edge_lengths).thisptr)))
            else:
                filename = tostring(edge_lengths)
                self.thisptr = new shared_ptr[Cpp_GillespieWorld](
                    new Cpp_GillespieWorld(filename))
        else:
            # XXX: GSLRandomNumberGenerator -> RandomNumberGenerator
            self.thisptr = new shared_ptr[Cpp_GillespieWorld](
                new Cpp_GillespieWorld(
                    deref((<Real3>edge_lengths).thisptr), deref(rng.thisptr)))

    def __dealloc__(self):
        # XXX: Here, we release shared pointer,
        #      and if reference count to the GillespieWorld object,
        #      it will be released automatically.
        del self.thisptr

    def set_t(self, Real t):
        """set_t(t)

        Set the value of the time of the world.

        Args:
            t (Real): the time of the world

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

    def volume(self):
        """Return the volume of the world."""
        return self.thisptr.get().volume()

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

    def list_particles(self, Species sp = None):
        """list_particles(sp) -> [(ParticleID, Particle)]

        Return the list of particles.
        A position of each particle is randomly generated.

        Args:
            sp (Species, optional): the species of particles to list up
                If no species is given, return the whole list of particles.

        Returns:
            list: the list of particles (of the given species)

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
        A position of each particle is randomly generated.

        Args:
            sp (Species): the species of particles to list up

        Returns:
            list: the list of particles of a given species

        """
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]] particles
        particles = self.thisptr.get().list_particles_exact(deref(sp.thisptr))

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

    def save(self, filename):
        """save(filename)

        Save self to a HDF5 file.

        Args:
            filename (str): a filename

        """
        self.thisptr.get().save(tostring(filename))

    def load(self, filename):
        """load(filename)

        Load self from a HDF5 file.

        Args:
            filename (str): a filename

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

cdef GillespieWorld GillespieWorld_from_Cpp_GillespieWorld(
    shared_ptr[Cpp_GillespieWorld] w):
    r = GillespieWorld(Real3(1, 1, 1))
    r.thisptr.swap(w)
    return r

## GillespieSimulator
#  a python wrapper for Cpp_GillespieSimulator
cdef class GillespieSimulator:
    """ A class running the simulation with the gillespie algorithm.

    GillespieSimulator(m, w)

    """

    def __init__(self, m, GillespieWorld w=None):
        """GillespieSimulator(m, w)
        GillespieSimulator(w)

        Constructor.

        Args:
            m (Model): A model
            w (GillespieWorld): A world

        """
        pass

    def __cinit__(self, m, GillespieWorld w=None):
        if w is None:
            self.thisptr = new Cpp_GillespieSimulator(
                deref((<GillespieWorld>m).thisptr))
        else:
            self.thisptr = new Cpp_GillespieSimulator(
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

    def set_t(self, Real new_t):
        """set_t(t)

        Set the current time.

        Args:
            t (Real): a current time.

        """
        self.thisptr.set_t(new_t)

    def set_dt(self, Real dt):
        """set_dt(dt)

        Set a step interval.

        Args:
            dt (Real): a step interval

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
        return GillespieWorld_from_Cpp_GillespieWorld(self.thisptr.world())

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

cdef GillespieSimulator GillespieSimulator_from_Cpp_GillespieSimulator(Cpp_GillespieSimulator* s):
    r = GillespieSimulator(
        Model_from_Cpp_Model(s.model()), GillespieWorld_from_Cpp_GillespieWorld(s.world()))
    del r.thisptr
    r.thisptr = s
    return r

## GillespieFactory
#  a python wrapper for Cpp_GillespieFactory
cdef class GillespieFactory:
    """ A factory class creating a GillespieWorld instance and a GillespieSimulator instance.

    GillespieFactory(GSLRandomNumberGenerator rng=None)

    """

    def __init__(self, GSLRandomNumberGenerator rng=None):
        """Constructor.

        Args:
            rng (GSLRandomNumberGenerator, optional): a random number generator.

        """
        pass

    def __cinit__(self, GSLRandomNumberGenerator rng=None):
        if rng is None:
            self.thisptr = new Cpp_GillespieFactory()
        else:
            self.thisptr = new Cpp_GillespieFactory(deref(rng.thisptr))

    def __dealloc__(self):
        del self.thisptr

    def create_world(self, arg1=None):
        """create_world(arg1=None) -> GillespieWorld

        Return a GillespieWorld instance.

        Args:
            arg1 (Real3): The lengths of edges of a GillespieWorld created

            or

            arg1 (str): The path of a HDF5 file for GillespieWorld

        Returns:
            GillespieWorld: the created world

        """
        if arg1 is None:
            return GillespieWorld_from_Cpp_GillespieWorld(
                shared_ptr[Cpp_GillespieWorld](self.thisptr.create_world()))
        elif isinstance(arg1, Real3):
            return GillespieWorld_from_Cpp_GillespieWorld(
                shared_ptr[Cpp_GillespieWorld](
                    self.thisptr.create_world(deref((<Real3>arg1).thisptr))))
        elif isinstance(arg1, str):
            return GillespieWorld_from_Cpp_GillespieWorld(
                shared_ptr[Cpp_GillespieWorld](self.thisptr.create_world(<string>(arg1))))
        else:
            return GillespieWorld_from_Cpp_GillespieWorld(
                shared_ptr[Cpp_GillespieWorld](self.thisptr.create_world(
                    deref(Cpp_Model_from_Model(arg1)))))

    def create_simulator(self, arg1, GillespieWorld arg2=None):
        """create_simulator(arg1, arg2) -> GillespieSimulator

        Return a GillespieSimulator instance.

        Args:
            arg1 (GillespieWorld): a world

            or

            arg1 (Model): a simulation model
            arg2 (GillespieWorld): a world

        Returns:
            GillespieSimulator: the created simulator

        """
        if arg2 is None:
            return GillespieSimulator_from_Cpp_GillespieSimulator(
                self.thisptr.create_simulator(deref((<GillespieWorld>arg1).thisptr)))
        else:
            return GillespieSimulator_from_Cpp_GillespieSimulator(
                self.thisptr.create_simulator(
                    deref(Cpp_Model_from_Model(arg1)), deref(arg2.thisptr)))

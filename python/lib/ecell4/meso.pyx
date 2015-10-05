import collections
from cython cimport address
from cython.operator cimport dereference as deref, preincrement as inc
from ecell4.core cimport *


## ReactionInfo
cdef class ReactionInfo:
    """A class stores detailed information about a reaction in meso.

    ReactionInfo(t, reactants, products, coord)

    """

    def __init__(self, Real t, reactants, products, coord):
        """Constructor.

        Args:
          t (Real): A time when a reaction occurred.
          reactants (list): A list of reactants.
            Reactants are given as a ``Species``.
          products (list): A list of products.
            Products are given as a ``Species``.
          coord (int): A coordinate where a reaction occurred.

        """
        pass  #XXX: only used for doc string


    def __cinit__(self, Real t, reactants, products, Integer coord):
        cdef vector[Cpp_Species] reactants_
        cdef vector[Cpp_Species] products_

        for sp in reactants:
            reactants_.push_back(deref((<Species>sp).thisptr))
        for sp in products:
            products_.push_back(deref((<Species>sp).thisptr))

        self.thisptr = new Cpp_ReactionInfo(t, reactants_, products_, coord)

    def __dealloc__(self):
        del self.thisptr

    def t(self):
        """Return a time when a reaction occurred."""
        return self.thisptr.t()

    def coordinate(self):
        """Return a coordinate where a reaction occurred."""
        return self.thisptr.coordinate()

    def reactants(self):
        """Return a list of reactants

        Returns:
            list: A list of ``Species``.

        """
        cdef vector[Cpp_Species] species = self.thisptr.reactants()

        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(
                 Species_from_Cpp_Species(
                     <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    def products(self):
        cdef vector[Cpp_Species] species = self.thisptr.products()
        """Return a list of products

        Returns:
            list: A list of ``Species``.

        """

        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(
                 Species_from_Cpp_Species(
                     <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

cdef ReactionInfo ReactionInfo_from_Cpp_ReactionInfo(Cpp_ReactionInfo* ri):
    cdef Cpp_ReactionInfo *new_obj = new Cpp_ReactionInfo(<Cpp_ReactionInfo> deref(ri))
    r = ReactionInfo(0, [], [], 0)
    del r.thisptr
    r.thisptr = new_obj
    return r

## MesoscopicWorld
#  a python wrapper for Cpp_MesoscopicWorld
cdef class MesoscopicWorld:

    def __cinit__(self, edge_lengths = None,
        matrix_sizes = None, GSLRandomNumberGenerator rng = None):
        cdef string filename

        if edge_lengths is None:
            self.thisptr = new shared_ptr[Cpp_MesoscopicWorld](new Cpp_MesoscopicWorld())
        elif matrix_sizes is None:
            if isinstance(edge_lengths, Real3):
                self.thisptr = new shared_ptr[Cpp_MesoscopicWorld](
                    new Cpp_MesoscopicWorld(deref((<Real3>edge_lengths).thisptr)))
            else:
                filename = tostring(edge_lengths)
                self.thisptr = new shared_ptr[Cpp_MesoscopicWorld](
                    new Cpp_MesoscopicWorld(filename))
        elif rng is None:
            if isinstance(matrix_sizes, Integer3):
                self.thisptr = new shared_ptr[Cpp_MesoscopicWorld](
                    new Cpp_MesoscopicWorld(
                        deref((<Real3>edge_lengths).thisptr),
                        deref((<Integer3>matrix_sizes).thisptr)))
            else:
                self.thisptr = new shared_ptr[Cpp_MesoscopicWorld](
                    new Cpp_MesoscopicWorld(
                        deref((<Real3>edge_lengths).thisptr),
                        <Real>matrix_sizes))
        else:
            if isinstance(matrix_sizes, Integer3):
                # XXX: GSLRandomNumberGenerator -> RandomNumberGenerator
                self.thisptr = new shared_ptr[Cpp_MesoscopicWorld](
                    new Cpp_MesoscopicWorld(
                        deref((<Real3>edge_lengths).thisptr),
                        deref((<Integer3>matrix_sizes).thisptr), deref(rng.thisptr)))
            else:
                self.thisptr = new shared_ptr[Cpp_MesoscopicWorld](
                    new Cpp_MesoscopicWorld(
                        deref((<Real3>edge_lengths).thisptr),
                        <Real>matrix_sizes, deref(rng.thisptr)))

    def __dealloc__(self):
        # XXX: Here, we release shared pointer,
        #      and if reference count to the MesoscopicWorld object,
        #      it will be released automatically.
        del self.thisptr

    def set_t(self, Real t):
        self.thisptr.get().set_t(t)

    def t(self):
        return self.thisptr.get().t()

    def edge_lengths(self):
        cdef Cpp_Real3 lengths = self.thisptr.get().edge_lengths()
        return Real3_from_Cpp_Real3(address(lengths))

    def matrix_sizes(self):
        cdef Cpp_Integer3 sizes = self.thisptr.get().matrix_sizes()
        return Integer3_from_Cpp_Integer3(address(sizes))

    def volume(self):
        return self.thisptr.get().volume()

    def subvolume(self):
        return self.thisptr.get().subvolume()

    def num_subvolumes(self, sp = None):
        if sp is None:
            return self.thisptr.get().num_subvolumes()
        else:
            return self.thisptr.get().num_subvolumes(deref((<Species>sp).thisptr))

    def subvolume_edge_lengths(self):
        cdef Cpp_Real3 lengths = self.thisptr.get().subvolume_edge_lengths()
        return Real3_from_Cpp_Real3(address(lengths))

    def global2coord(self, Integer3 g):
        return self.thisptr.get().global2coord(deref(g.thisptr))

    def coord2global(self, Integer c):
        cdef Cpp_Integer3 g = self.thisptr.get().coord2global(c)
        return Integer3_from_Cpp_Integer3(address(g))

    def position2global(self, Real3 pos):
        cdef Cpp_Integer3 g = self.thisptr.get().position2global(deref(pos.thisptr))
        return Integer3_from_Cpp_Integer3(address(g))

    def position2coordinate(self, Real3 pos):
        return self.thisptr.get().position2coordinate(deref(pos.thisptr))

    def num_molecules(self, Species sp, c = None):
        if c is None:
            return self.thisptr.get().num_molecules(deref(sp.thisptr))
        elif isinstance(c, Integer3):
            return self.thisptr.get().num_molecules(deref(sp.thisptr), deref((<Integer3>c).thisptr))
        else:
            return self.thisptr.get().num_molecules(deref(sp.thisptr), <Integer>c)

    def num_molecules_exact(self, Species sp, c = None):
        if c is None:
            return self.thisptr.get().num_molecules_exact(deref(sp.thisptr))
        elif isinstance(c, Integer3):
            return self.thisptr.get().num_molecules_exact(deref(sp.thisptr), deref((<Integer3>c).thisptr))
        else:
            return self.thisptr.get().num_molecules_exact(deref(sp.thisptr), <Integer>c)

    def add_molecules(self, Species sp, Integer num, c = None):
        if c is None:
            self.thisptr.get().add_molecules(deref(sp.thisptr), num)
        elif isinstance(c, Integer3):
            self.thisptr.get().add_molecules(deref(sp.thisptr), num, deref((<Integer3>c).thisptr))
        elif hasattr(c, "as_base"):
            self.thisptr.get().add_molecules(
                deref(sp.thisptr), num, deref((<Shape>(c.as_base())).thisptr))
        else:
            self.thisptr.get().add_molecules(deref(sp.thisptr), num, <Integer>c)

    def remove_molecules(self, Species sp, Integer num, c = None):
        if c is None:
            self.thisptr.get().remove_molecules(deref(sp.thisptr), num)
        elif isinstance(c, Integer3):
            self.thisptr.get().remove_molecules(deref(sp.thisptr), num, deref((<Integer3>c).thisptr))
        else:
            self.thisptr.get().remove_molecules(deref(sp.thisptr), num, <Integer>c)

    def add_structure(self, Species sp, shape):
        self.thisptr.get().add_structure(
            deref(sp.thisptr), deref((<Shape>(shape.as_base())).thisptr))

    def get_volume(self, Species sp):
        return self.thisptr.get().get_volume(deref(sp.thisptr))

    def has_structure(self, Species sp):
        return self.thisptr.get().has_structure(deref(sp.thisptr))

    def on_structure(self, Species sp, Integer3 g):
        return self.thisptr.get().on_structure(deref(sp.thisptr), deref(g.thisptr))

    def check_structure(self, Species sp, Integer3 g):
        return self.thisptr.get().check_structure(deref(sp.thisptr), deref(g.thisptr))

    def get_occupancy(self, Species sp, g):
        if isinstance(g, Integer3):
            return self.thisptr.get().get_occupancy(
                deref(sp.thisptr), deref((<Integer3>g).thisptr))
        else:
            return self.thisptr.get().get_occupancy(deref(sp.thisptr), <Integer>g)

    def list_species(self):
        cdef vector[Cpp_Species] species = self.thisptr.get().list_species()

        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(
                 Species_from_Cpp_Species(
                     <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    def list_coordinates(self, Species sp):
        return self.thisptr.get().list_coordinates(deref(sp.thisptr))

    def list_coordinates_exact(self, Species sp):
        return self.thisptr.get().list_coordinates_exact(deref(sp.thisptr))

    def list_particles(self, Species sp = None):
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
        self.thisptr.get().save(tostring(filename))

    def load(self, filename):
        self.thisptr.get().load(tostring(filename))

    def bind_to(self, m):
        self.thisptr.get().bind_to(deref(Cpp_Model_from_Model(m)))

    def rng(self):
        return GSLRandomNumberGenerator_from_Cpp_RandomNumberGenerator(
            self.thisptr.get().rng())

    def as_base(self):
        retval = Space()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Space](
            <shared_ptr[Cpp_Space]>deref(self.thisptr))
        return retval

cdef MesoscopicWorld MesoscopicWorld_from_Cpp_MesoscopicWorld(
    shared_ptr[Cpp_MesoscopicWorld] w):
    r = MesoscopicWorld(Real3(1, 1, 1), Integer3(1, 1, 1))
    r.thisptr.swap(w)
    return r

## MesoscopicSimulator
#  a python wrapper for Cpp_MesoscopicSimulator
cdef class MesoscopicSimulator:

    def __cinit__(self, m, MesoscopicWorld w=None):
        if w is None:
            self.thisptr = new Cpp_MesoscopicSimulator(
                deref((<MesoscopicWorld>m).thisptr))
        else:
            self.thisptr = new Cpp_MesoscopicSimulator(
                deref(Cpp_Model_from_Model(m)), deref(w.thisptr))

    def __dealloc__(self):
        del self.thisptr

    def num_steps(self):
        return self.thisptr.num_steps()

    def step(self, upto = None):
        if upto is None:
            self.thisptr.step()
        else:
            return self.thisptr.step(<Real> upto)

    def t(self):
        return self.thisptr.t()

    def dt(self):
        return self.thisptr.dt()

    def next_time(self):
        return self.thisptr.next_time()

    def check_reaction(self):
        """Return if any reaction occurred at the last step, or not."""
        return self.thisptr.check_reaction()

    def last_reactions(self):
        """Return a list of reactions, which occurred at the last step.

        Returns:
            list: A list of pairs of ``ReactionRule`` and ``ReactionInfo``.

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
        self.thisptr.set_t(new_t)

    def set_dt(self, Real dt):
        self.thisptr.set_dt(dt)

    def initialize(self):
        self.thisptr.initialize()

    def model(self):
        return Model_from_Cpp_Model(self.thisptr.model())

    def world(self):
        return MesoscopicWorld_from_Cpp_MesoscopicWorld(self.thisptr.world())

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

cdef MesoscopicSimulator MesoscopicSimulator_from_Cpp_MesoscopicSimulator(
    Cpp_MesoscopicSimulator* s):
    r = MesoscopicSimulator(
        Model_from_Cpp_Model(s.model()),
        MesoscopicWorld_from_Cpp_MesoscopicWorld(s.world()))
    del r.thisptr
    r.thisptr = s
    return r

## MesoscopicFactory
#  a python wrapper for Cpp_MesoscopicFactory
cdef class MesoscopicFactory:

    def __cinit__(self, matrix_sizes=None, GSLRandomNumberGenerator rng=None):
        if rng is not None:
            if isinstance(matrix_sizes, Integer3):
                self.thisptr = new Cpp_MesoscopicFactory(
                    deref((<Integer3>matrix_sizes).thisptr), deref(rng.thisptr))
            else:
                self.thisptr = new Cpp_MesoscopicFactory(<Real>matrix_sizes, deref(rng.thisptr))
        elif matrix_sizes is not None:
            if isinstance(matrix_sizes, Integer3):
                self.thisptr = new Cpp_MesoscopicFactory(
                    deref((<Integer3>matrix_sizes).thisptr))
            else:
                self.thisptr = new Cpp_MesoscopicFactory(<Real>matrix_sizes)
        else:
            self.thisptr = new Cpp_MesoscopicFactory()

    def __dealloc__(self):
        del self.thisptr

    def create_world(self, arg1):
        if isinstance(arg1, Real3):
            return MesoscopicWorld_from_Cpp_MesoscopicWorld(
                shared_ptr[Cpp_MesoscopicWorld](
                    self.thisptr.create_world(deref((<Real3>arg1).thisptr))))
        elif isinstance(arg1, str):
            return MesoscopicWorld_from_Cpp_MesoscopicWorld(
                shared_ptr[Cpp_MesoscopicWorld](self.thisptr.create_world(<string>(arg1))))
        else:
            return MesoscopicWorld_from_Cpp_MesoscopicWorld(
                shared_ptr[Cpp_MesoscopicWorld](self.thisptr.create_world(
                    deref(Cpp_Model_from_Model(arg1)))))

    def create_simulator(self, arg1, MesoscopicWorld arg2=None):
        if arg2 is None:
            return MesoscopicSimulator_from_Cpp_MesoscopicSimulator(
                self.thisptr.create_simulator(deref((<MesoscopicWorld>arg1).thisptr)))
        else:
            return MesoscopicSimulator_from_Cpp_MesoscopicSimulator(
                self.thisptr.create_simulator(
                    deref(Cpp_Model_from_Model(arg1)), deref(arg2.thisptr)))

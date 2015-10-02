import collections
from cython cimport address
from cython.operator cimport dereference as deref, preincrement as inc
from ecell4.core cimport *


## ReactionInfo
cdef class ReactionInfo:
    """A class stores detailed information about a reaction in gillespie.

    ReactionInfo(t, reactants, products)

    """

    def __init__(self, Real t, reactants, products):
        """Constructor.

        Args:
          t (Real): A time when a reaction occurred
          reactants (list): A list of reactants.
            Reactants are given as a ``Species``.
          products (list): A list of products.
            Products are given as a ``Species``.

        """
        pass  #XXX: only used for doc string

    def __cinit__(self, Real t, reactants, products):
        cdef vector[Cpp_Species] reactants_
        cdef vector[Cpp_Species] products_

        for sp in reactants:
            reactants_.push_back(deref((<Species>sp).thisptr))
        for sp in products:
            products_.push_back(deref((<Species>sp).thisptr))

        self.thisptr = new Cpp_ReactionInfo(t, reactants_, products_)

    def __dealloc__(self):
        del self.thisptr

    def t(self):
        """Return a time when a reaction occurred."""
        return self.thisptr.t()

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
        """Return a list of products

        Returns:
            list: A list of ``Species``.

        """
        cdef vector[Cpp_Species] species = self.thisptr.products()

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
    r = ReactionInfo(0, [], [])
    del r.thisptr
    r.thisptr = new_obj
    return r

## GillespieWorld
#  a python wrapper for Cpp_GillespieWorld
cdef class GillespieWorld:

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

    def add_molecules(self, Species sp, Integer num, shape=None):
        if shape is None:
            self.thisptr.get().add_molecules(deref(sp.thisptr), num)
        else:
            self.thisptr.get().add_molecules(
                deref(sp.thisptr), num, deref((<Shape>(shape.as_base())).thisptr))

    def remove_molecules(self, Species sp, Integer num):
        self.thisptr.get().remove_molecules(deref(sp.thisptr), num)

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

cdef GillespieWorld GillespieWorld_from_Cpp_GillespieWorld(
    shared_ptr[Cpp_GillespieWorld] w):
    r = GillespieWorld(Real3(1, 1, 1))
    r.thisptr.swap(w)
    return r

## GillespieSimulator
#  a python wrapper for Cpp_GillespieSimulator
cdef class GillespieSimulator:

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
        return GillespieWorld_from_Cpp_GillespieWorld(self.thisptr.world())

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

cdef GillespieSimulator GillespieSimulator_from_Cpp_GillespieSimulator(Cpp_GillespieSimulator* s):
    r = GillespieSimulator(
        Model_from_Cpp_Model(s.model()), GillespieWorld_from_Cpp_GillespieWorld(s.world()))
    del r.thisptr
    r.thisptr = s
    return r

## GillespieFactory
#  a python wrapper for Cpp_GillespieFactory
cdef class GillespieFactory:

    def __cinit__(self, GSLRandomNumberGenerator rng=None):
        if rng is None:
            self.thisptr = new Cpp_GillespieFactory()
        else:
            self.thisptr = new Cpp_GillespieFactory(deref(rng.thisptr))

    def __dealloc__(self):
        del self.thisptr

    def create_world(self, arg1):
        if isinstance(arg1, Real3):
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
        if arg2 is None:
            return GillespieSimulator_from_Cpp_GillespieSimulator(
                self.thisptr.create_simulator(deref((<GillespieWorld>arg1).thisptr)))
        else:
            return GillespieSimulator_from_Cpp_GillespieSimulator(
                self.thisptr.create_simulator(
                    deref(Cpp_Model_from_Model(arg1)), deref(arg2.thisptr)))

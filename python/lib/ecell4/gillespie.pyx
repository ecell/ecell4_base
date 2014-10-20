from cython cimport address
from cython.operator cimport dereference as deref, preincrement as inc


## GillespieWorld
#  a python wrapper for Cpp_GillespieWorld
cdef class GillespieWorld:

    def __cinit__(self, edge_lengths = None, GSLRandomNumberGenerator rng = None):
        cdef string filename

        if edge_lengths is None:
            self.thisptr = new shared_ptr[Cpp_GillespieWorld](new Cpp_GillespieWorld())
        elif rng is None:
            if isinstance(edge_lengths, Position3):
                self.thisptr = new shared_ptr[Cpp_GillespieWorld](
                    new Cpp_GillespieWorld(deref((<Position3>edge_lengths).thisptr)))
            else:
                filename = edge_lengths
                self.thisptr = new shared_ptr[Cpp_GillespieWorld](
                    new Cpp_GillespieWorld(filename))
        else:
            # XXX: GSLRandomNumberGenerator -> RandomNumberGenerator
            self.thisptr = new shared_ptr[Cpp_GillespieWorld](
                new Cpp_GillespieWorld(
                    deref((<Position3>edge_lengths).thisptr), deref(rng.thisptr)))

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
        cdef Cpp_Position3 lengths = self.thisptr.get().edge_lengths()
        return Position3_from_Cpp_Position3(address(lengths))

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

    def save(self, string filename):
        self.thisptr.get().save(filename)

    def load(self, string filename):
        self.thisptr.get().load(filename)

    def bind_to(self, m):
        self.thisptr.get().bind_to(deref(Cpp_Model_from_Model(m)))

    def rng(self):
        return GSLRandomNumberGenerator_from_Cpp_RandomNumberGenerator(
            self.thisptr.get().rng())

cdef GillespieWorld GillespieWorld_from_Cpp_GillespieWorld(
    shared_ptr[Cpp_GillespieWorld] w):
    r = GillespieWorld(Position3(1, 1, 1))
    r.thisptr.swap(w)
    return r

## GillespieSimulator
#  a python wrapper for Cpp_GillespieSimulator
cdef class GillespieSimulator:

    def __cinit__(self, m, GillespieWorld w):
        if w is None:
            self.thisptr = new Cpp_GillespieSimulator(
                deref((<GillespieWorld>w).thisptr))
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

    def last_reactions(self):
        cdef vector[Cpp_ReactionRule] reactions = self.thisptr.last_reactions()
        cdef vector[Cpp_ReactionRule].iterator it = reactions.begin()
        retval = []
        while it != reactions.end():
            retval.append(ReactionRule_from_Cpp_ReactionRule(
                <Cpp_ReactionRule*>(address(deref(it)))))
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
        else:
            for obs in observers:
                tmp.push_back(deref((<Observer>(obs.as_base())).thisptr))
            self.thisptr.run(duration, tmp)

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
        if isinstance(arg1, Position3):
            return GillespieWorld_from_Cpp_GillespieWorld(
                shared_ptr[Cpp_GillespieWorld](
                    self.thisptr.create_world(deref((<Position3>arg1).thisptr))))
        else:
            return GillespieWorld_from_Cpp_GillespieWorld(
                shared_ptr[Cpp_GillespieWorld](self.thisptr.create_world(<string>(arg1))))

    def create_simulator(self, arg1, GillespieWorld arg2=None):
        if arg2 is None:
            return GillespieSimulator_from_Cpp_GillespieSimulator(
                self.thisptr.create_simulator(deref((<GillespieWorld>arg1).thisptr)))
        else:
            return GillespieSimulator_from_Cpp_GillespieSimulator(
                self.thisptr.create_simulator(
                    deref(Cpp_Model_from_Model(arg1)), deref(arg2.thisptr)))

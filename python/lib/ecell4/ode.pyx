from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address
from libcpp.string cimport string
from libcpp.vector cimport vector

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *

## ODEWorld
#  a python wrapper for Cpp_ODEWorld
cdef class ODEWorld:

    def __cinit__(self, Position3 edge_lengths):
        # XXX: GSLRandomNumberGenerator -> RandomNumberGenerator
        self.thisptr = new shared_ptr[Cpp_ODEWorld](
            new Cpp_ODEWorld(deref(edge_lengths.thisptr)))

    def __dealloc__(self):
        # XXX: Here, we release shared pointer,
        #      and if reference count to the ODEWorld object become zero,
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

    def list_species(self):
        cdef vector[Cpp_Species] raw_list_species = self.thisptr.get().list_species()
        retval = []
        cdef vector[Cpp_Species].iterator it = raw_list_species.begin()
        while it != raw_list_species.end():
            retval.append(
                Species_from_Cpp_Species(<Cpp_Species*> (address(deref(it)))))
            inc(it)
        return retval

    def set_volume(self, Real vol):
        self.thisptr.get().set_volume(vol)

    def add_molecules(self, Species sp, Real num):
        self.thisptr.get().add_molecules(deref(sp.thisptr), num)

    def remove_molecules(self, Species sp, Real num):
        self.thisptr.get().remove_molecules(deref(sp.thisptr), num)

    def get_value(self, Species sp):
        return self.thisptr.get().get_value(deref(sp.thisptr))

    def set_value(self, Species sp, Real num):
        self.thisptr.get().set_value(deref(sp.thisptr), num)

    def save(self, string filename):
        self.thisptr.get().save(filename)

    def load(self, string filename):
        self.thisptr.get().load(filename)

    def has_species(self, Species sp):
        return self.thisptr.get().has_species(deref(sp.thisptr))

    def reserve_species(self, Species sp):
        self.thisptr.get().reserve_species(deref(sp.thisptr))

    def release_species(self, Species sp):
        self.thisptr.get().release_species(deref(sp.thisptr))

    def bind_to(self, NetworkModel m):
        self.thisptr.get().bind_to(deref(m.thisptr))

cdef ODEWorld ODEWorld_from_Cpp_ODEWorld(
    shared_ptr[Cpp_ODEWorld] w):
    r = ODEWorld(Position3(1, 1, 1))
    r.thisptr.swap(w)
    return r

## ODESimulator
#  a python wrapper for Cpp_ODESimulator
cdef class ODESimulator:

    def __cinit__(self, NetworkModel m, ODEWorld w):
        self.thisptr = new Cpp_ODESimulator(
            deref(m.thisptr), deref(w.thisptr))

    def __dealloc__(self):
        del self.thisptr

    def num_steps(self):
        return self.thisptr.num_steps()

    def next_time(self):
        return self.thisptr.next_time()

    def dt(self):
        return self.thisptr.dt()

    def step(self, upto = None):
        if upto is None:
            self.thisptr.step()
        else:
            return self.thisptr.step(upto)

    def t(self):
        return self.thisptr.t()

    def set_t(self, Real new_t):
        self.thisptr.set_t(new_t)

    def set_dt(self, Real dt):
        self.thisptr.set_dt(dt)

    def initialize(self):
        self.thisptr.initialize()

    def last_reactions(self):
        cdef vector[Cpp_ReactionRule] reactions = self.thisptr.last_reactions()
        cdef vector[Cpp_ReactionRule].iterator it = reactions.begin()
        retval = []
        while it != reactions.end():
            retval.append(ReactionRule_from_Cpp_ReactionRule(
                <Cpp_ReactionRule*>(address(deref(it)))))
            inc(it)
        return retval

    def model(self):
        return NetworkModel_from_Cpp_NetworkModel(self.thisptr.model())

    def world(self):
        return ODEWorld_from_Cpp_ODEWorld(self.thisptr.world())

    def run(self, Real duration, observers=None):
        cdef vector[shared_ptr[Cpp_Observer]] tmp

        if observers is None:
            self.thisptr.run(duration)
        else:
            for obs in observers:
                tmp.push_back(deref((<Observer>(obs.as_base())).thisptr))
            self.thisptr.run(duration, tmp)

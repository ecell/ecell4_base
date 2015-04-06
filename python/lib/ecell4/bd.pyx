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
        cdef pair[pair[Cpp_ParticleID, Cpp_Particle], bool] retval

        if arg2 is None:
            retval = self.thisptr.get().new_particle(deref((<Particle> arg1).thisptr))
        else:
            retval = self.thisptr.get().new_particle(deref((<Species> arg1).thisptr), deref(arg2.thisptr))
        return ((ParticleID_from_Cpp_ParticleID(address(retval.first.first)), Particle_from_Cpp_Particle(address(retval.first.second))), retval.second)

    def set_t(self, Real t):
        self.thisptr.get().set_t(t)

    def t(self):
        return self.thisptr.get().t()

    def edge_lengths(self):
        cdef Cpp_Real3 lengths = self.thisptr.get().edge_lengths()
        return Real3_from_Cpp_Real3(address(lengths))

    def num_particles(self, Species sp = None):
        if sp is None:
            return self.thisptr.get().num_particles()
        else:
            return self.thisptr.get().num_particles(deref(sp.thisptr))

    def num_particles_exact(self, Species sp):
        return self.thisptr.get().num_particles_exact(deref(sp.thisptr))

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
        return self.thisptr.get().has_particle(deref(pid.thisptr))

    def update_particle(self, ParticleID pid, Particle p):
        return self.thisptr.get().update_particle(deref(pid.thisptr), deref(p.thisptr))

    def get_particle(self, ParticleID pid):
        cdef pair[Cpp_ParticleID, Cpp_Particle] \
            pid_particle_pair = self.thisptr.get().get_particle(deref(pid.thisptr))
        return (ParticleID_from_Cpp_ParticleID(address(pid_particle_pair.first)),
                Particle_from_Cpp_Particle(address(pid_particle_pair.second)))

    def remove_particle(self, ParticleID pid):
        self.thisptr.get().remove_particle(deref(pid.thisptr))

    def list_particles_within_radius(
        self, Real3 pos, Real radius,
        ParticleID ignore1 = None, ParticleID ignore2 = None):
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
        cdef Cpp_Real3 newpos = self.thisptr.get().periodic_transpose(
            deref(pos1.thisptr), deref(pos2.thisptr))
        return Real3_from_Cpp_Real3(address(newpos))

    def apply_boundary(self, Real3 pos):
        cdef Cpp_Real3 newpos = self.thisptr.get().apply_boundary(deref(pos.thisptr))
        return Real3_from_Cpp_Real3(address(newpos))

    def distance_sq(self, Real3 pos1, Real3 pos2):
        return self.thisptr.get().distance_sq(deref(pos1.thisptr), deref(pos2.thisptr))

    def distance(self, Real3 pos1, Real3 pos2):
        return self.thisptr.get().distance(deref(pos1.thisptr), deref(pos2.thisptr))

    def volume(self):
        return self.thisptr.get().volume()

    # def has_species(self, Species sp):
    #     return self.thisptr.get().has_species(deref(sp.thisptr))

    def num_molecules(self, Species sp):
        return self.thisptr.get().num_molecules(deref(sp.thisptr))

    def num_molecules_exact(self, Species sp):
        return self.thisptr.get().num_molecules_exact(deref(sp.thisptr))

    # def add_species(self, Species sp):
    #     self.thisptr.get().add_species(deref(sp.thisptr))

    # def add_molecules(self, Species sp, Integer num):
    #     self.thisptr.get().add_molecules(deref(sp.thisptr), num)

    def add_molecules(self, Species sp, Integer num, shape=None):
        if shape is None:
            self.thisptr.get().add_molecules(deref(sp.thisptr), num)
        else:
            self.thisptr.get().add_molecules(
                deref(sp.thisptr), num, deref((<Shape>(shape.as_base())).thisptr))

    def remove_molecules(self, Species sp, Integer num):
        self.thisptr.get().remove_molecules(deref(sp.thisptr), num)

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

cdef BDWorld BDWorld_from_Cpp_BDWorld(
    shared_ptr[Cpp_BDWorld] w):
    r = BDWorld(Real3(1, 1, 1))
    r.thisptr.swap(w)
    return r

## BDSimulator
#  a python wrapper for Cpp_BDSimulator
cdef class BDSimulator:

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
        return self.thisptr.num_steps()

    def step(self, upto = None):
        if upto is None:
            self.thisptr.step()
        else:
            return self.thisptr.step(upto)

    def t(self):
        return self.thisptr.t()

    def dt(self):
        return self.thisptr.dt()

    def set_dt(self, Real& dt):
        self.thisptr.set_dt(dt)

    def next_time(self):
        return self.thisptr.next_time()

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
        return Model_from_Cpp_Model(self.thisptr.model())

    def world(self):
        return BDWorld_from_Cpp_BDWorld(self.thisptr.world())

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

cdef BDSimulator BDSimulator_from_Cpp_BDSimulator(Cpp_BDSimulator* s):
    r = BDSimulator(
        Model_from_Cpp_Model(s.model()), BDWorld_from_Cpp_BDWorld(s.world()))
    del r.thisptr
    r.thisptr = s
    return r

## BDFactory
#  a python wrapper for Cpp_BDFactory
cdef class BDFactory:

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
        if isinstance(arg1, Real3):
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
        if arg2 is None:
            return BDSimulator_from_Cpp_BDSimulator(
                self.thisptr.create_simulator(deref((<BDWorld>arg1).thisptr)))
        else:
            return BDSimulator_from_Cpp_BDSimulator(
                self.thisptr.create_simulator(
                    deref(Cpp_Model_from_Model(arg1)), deref(arg2.thisptr)))

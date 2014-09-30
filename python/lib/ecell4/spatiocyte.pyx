from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address
from libcpp.string cimport string
from libcpp.vector cimport vector

from ecell4.types cimport *
from ecell4.shared_ptr cimport shared_ptr
from ecell4.core cimport *


## SpatiocyteWorld
#  a python wrapper for Cpp_SpatiocyteWorld
cdef class SpatiocyteWorld:

    def __cinit__(self, Position3 edge_lengths, Real voxel_radius):
        self.thisptr = new shared_ptr[Cpp_SpatiocyteWorld](
            new Cpp_SpatiocyteWorld(deref(edge_lengths.thisptr), voxel_radius))

    def __dealloc__(self):
        # XXX: Here, we release shared pointer,
        #      and if reference count to the SpatiocyteWorld object,
        #      it will be released automatically.
        del self.thisptr

    def set_t(self, Real t):
        self.thisptr.get().set_t(t)

    def t(self):
        return self.thisptr.get().t()

    def dt(self):
        return self.thisptr.get().dt()

    def voxel_radius(self):
        return self.thisptr.get().voxel_radius()

    def num_voxels(self):
        return self.thisptr.get().num_voxels()

    def get_voxel(self, Integer coord):
        cdef Cpp_Voxel voxel = self.thisptr.get().get_voxel(coord)
        return Voxel_from_Cpp_Voxel(address(voxel))

    def coordinates(self, Species sp):
        return self.thisptr.get().coordinates(deref(sp.thisptr))

    def edge_lengths(self):
        cdef Cpp_Position3 lengths = self.thisptr.get().edge_lengths()
        return Position3_from_Cpp_Position3(address(lengths))

    def num_particles(self, Species sp):
        if sp is None:
            return self.thisptr.get().num_particles()
        else:
            return self.thisptr.get().num_particles(deref(sp.thisptr))

    def list_particles(self, Species sp):
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]] \
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

    def volume(self):
        return self.thisptr.get().volume()

    # def num_species(self):
    #     return self.thisptr.get().num_species()

    def has_species(self, Species sp):
        return self.thisptr.get().has_species(deref(sp.thisptr))

    def num_molecules(self, Species sp):
        return self.thisptr.get().num_molecules(deref(sp.thisptr))

    # def add_species(self, Species sp):
    #     self.thisptr.get().add_species(deref(sp.thisptr))

    def add_molecules(self, Species sp, Integer num):
        self.thisptr.get().add_molecules(deref(sp.thisptr), num)

    def save(self, string filename):
        self.thisptr.get().save(filename)

## SpatiocyteSimulator
#  a python wrapper for Cpp_SpatiocyteSimulator
cdef class SpatiocyteSimulator:

    def __cinit__(self, NetworkModel m, SpatiocyteWorld w):
        self.thisptr = new Cpp_SpatiocyteSimulator(
            deref(m.thisptr), deref(w.thisptr))

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

    def initialize(self):
        self.thisptr.initialize()

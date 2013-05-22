from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address

from libcpp.pair cimport pair
from libcpp.vector cimport vector


cdef class ParticleSpaceVectorImpl:

    def __cinit__(self, Position3 edge_lengths):
        self.thisptr = new Cpp_ParticleSpaceVectorImpl(deref(edge_lengths.thisptr))

    def __dealloc__(self):
        del self.thisptr

    def edge_lengths(self):
        cdef Cpp_Position3 r = self.thisptr.edge_lengths()
        return Position3_from_Cpp_Position3(address(r))

    def update_particle(self, ParticleID pid, Particle p):
        return self.thisptr.update_particle(deref(pid.thisptr), deref(p.thisptr))

    def list_particles(self):
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]] \
            particles = self.thisptr.list_particles()
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

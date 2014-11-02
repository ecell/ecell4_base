from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address

from libcpp.pair cimport pair
from libcpp.vector cimport vector


cdef class ParticleSpaceVectorImpl:

    def __cinit__(self, Real3 edge_lengths):
        self.thisptr = new Cpp_ParticleSpaceVectorImpl(deref(edge_lengths.thisptr))

    def __dealloc__(self):
        del self.thisptr

    def edge_lengths(self):
        cdef Cpp_Real3 r = self.thisptr.edge_lengths()
        return Real3_from_Cpp_Real3(address(r))

    def num_particles(self, Species sp = None):
        if sp is None:
            return self.thisptr.num_particles()
        else:
            return self.thisptr.num_particles(deref(sp.thisptr))

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

    def has_particle(self, ParticleID pid):
        return self.thisptr.has_particle(deref(pid.thisptr))

    def get_particle(self, ParticleID pid):
        cdef pair[Cpp_ParticleID, Cpp_Particle] raw_particle_set = self.thisptr.get_particle(deref(pid.thisptr))
        return ( ParticleID_from_Cpp_ParticleID(address(raw_particle_set.first)), Particle_from_Cpp_Particle(address(raw_particle_set.second)) )

    def remove_particle(self, ParticleID pid):
        self.thisptr.remove_particle(deref(pid.thisptr))

    def list_particles_within_radius(self, Real3 pos, Real radius, ParticleID ignore1 = None, ParticleID ignore2 = None):
        cdef vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real] ] raw_list_particles_within_radius
        if ignore1 == None and ignore2 == None:
            raw_list_particles_within_radius = self.thisptr.list_particles_within_radius(
                    deref(pos.thisptr), radius)
        elif ignore2 == None:
            raw_list_particles_within_radius = self.thisptr.list_particles_within_radius(
                    deref(pos.thisptr), radius, deref(ignore1.thisptr) )
        else:
            raw_list_particles_within_radius = self.thisptr.list_particles_within_radius(
                    deref(pos.thisptr), radius, deref(ignore1.thisptr), deref(ignore2.thisptr))

        retval = []
        cdef vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real] ].iterator it = raw_list_particles_within_radius.begin()
        while it != raw_list_particles_within_radius.end():
            retval.append( 
                ((ParticleID_from_Cpp_ParticleID(
                     <Cpp_ParticleID*>(address(deref(it).first.first))),
                 Particle_from_Cpp_Particle(
                     <Cpp_Particle*>(address(deref(it).first.second)))),
                 deref(it).second)
                )
        return retval


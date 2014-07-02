# coding: utf-8
"""ec4vis.plugins.particle_space --- Draft implementation of ParticleSpace.
"""


class Particle(object):

    def __init__(self, sid, pos, radius, D = 0.0):
        self.sid = sid
        self.position = pos
        self.radius = radius
        self.D = D

# end of Particle

class ParticleSpace(object):

    def __init__(self):
        self.__species_pool = {}

    @property
    def species(self):
        if self.__species_pool is None:
            return []
        else:
            return self.__species_pool.keys()

    def add_particle(self, pid, particle):
        if particle.sid not in self.__species_pool.keys():
            self.__species_pool[particle.sid] = [(pid, particle)]
        else:
            self.__species_pool[particle.sid].append((pid, particle))

    def list_particles(self, sid=None):
        if sid is None:
            retval = []
            for sid, particles in self.__species_pool.items():
                retval.extend(particles)
            return retval
        elif sid not in self.__species_pool:
            return None
        return self.__species_pool[sid]

    def num_particles(self, sid=None):
        if sid is None:
            counts = [len(particles) for particles in self.__species_pool.values()]
            return sum(counts)
        if sid not in self.__species_pool:
            return 0
        return len(self.__species_pool[sid])

# end of ParticleSpace

# coding: utf-8
"""ec4vis.plugins.particle_space --- Draft implementation of ParticleSpace.
"""

import numpy

# try:
#     import ec4vis
# except ImportError:
#     import sys, os
#     p = os.path.abspath(__file__); sys.path.insert(0, p[: p.rindex(os.sep + 'ec4vis')])

from spatiocyte_tools import coord2point
from particle_space import Particle

class LatticeParticle(object):

    def __init__(self, sid, coord):
        self.sid = sid
        self.coord = coord

# end of LatticeParticle

class LatticeParticleSpace(object):

    def __init__(self, col_size, row_size, layer_size, lspecies, voxel_radius):
        self.__lattice_pool = {}
        #self.__col_size = col_size
        self.__row_size = row_size
        self.__layer_size = layer_size
        self.__lspecies = lspecies
        self.__voxel_radius = voxel_radius

    @property
    def species(self):
        if self.__lattice_pool is None:
            return []
        else:
            species = []
            for key in self.__lattice_pool.keys():
                (string, radius) = self.__lspecies[key]
                print "key : %d, string : %s, radius : %f" % (key, string, radius)
                species.append(string)
            return species

    def add_particle(self, lattice):
        sid = lattice.sid
        if sid not in self.__lattice_pool.keys():
            self.__lattice_pool[sid] = [lattice]
        else:
            self.__lattice_pool[sid].append(lattice)

    def __lattice2particle(self, lattice):
        pos = coord2point(lattice.coord, self.__row_size, self.__layer_size)
        pos = numpy.array(pos) * 2 * self.__voxel_radius
        (string, radius) = self.__lspecies[lattice.sid]
        return Particle(string, pos, radius)

    def __lattices2particles(self, lattices):
        particles = []
        pid = 0
        for lattice in lattices:
            particles.append((pid,self.__lattice2particle(lattice)))
        return particles

    def __string2key(self, string):
        for key in range(len(self.__lspecies)):
            (s, r) = self.__lspecies[key]
            if (string is s):
                return key
        return None

    def list_particles(self, sid=None):
        if sid is None:
            retval = []
            for lattices in self.__lattice_pool.values():
                retval.extend(self.__lattices2particles(lattices))
            return retval
        else:
            sid = self.__string2key(sid)
            if sid is None:
                return None
            return self.__lattices2particles(self.__lattice_pool[sid])

    def num_particles(self, sid=None):
        if sid is None:
            counts = [len(lattices) for lattices in self.__lattice_pool.values()]
            return sum(counts)
        else:
            sid = self.__string2key(sid)
            if sid is None:
                return 0
            return len(self.__lattice_pool[sid])

# end of LatticeParticleSpace

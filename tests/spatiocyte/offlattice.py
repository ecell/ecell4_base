import unittest
from ecell4_base.core import *
from ecell4_base.spatiocyte import *

class OffLatticeTest(unittest.TestCase):

    def setUp(self):
        self.voxel_radius = 0.005

        coordinates = [Real3(x, 0.0, 0.0) for x in range(0,10)]
        connections = [(x, x+1) for x in range(0,9)]
        self.offlattice = OffLattice(self.voxel_radius, coordinates, connections)

    def test_constructor(self):
        species = Species('Base')
        world = SpatiocyteWorld(ones(), self.voxel_radius)
        world.add_offlattice(species, self.offlattice)

import unittest
from ecell4_base.core import *
from ecell4_base.spatiocyte import *

class OffLatticeTest(unittest.TestCase):

    def setUp(self):
        self.coordinates = [Real3(x, 0.0, 0.0) for x in range(0,10)]
        self.connections = [(x, x+1) for x in range(0,9)]

    def test_constructor(self):
        voxel_radius = 0.005
        offlattice = OffLatticeSpace(voxel_radius, self.coordinates, self.connections)
        world = SpatiocyteWorld(ones(), voxel_radius)
        world.add_space(offlattice)

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

    def test_add_molecules(self):
        base = Species('Base', radius=self.voxel_radius, D=0.0)
        species = Species('A', radius=self.voxel_radius, D=1.0, location='Base')

        model = NetworkModel()
        model.add_species_attribute(base)
        model.add_species_attribute(species)

        world = SpatiocyteWorld(ones(), self.voxel_radius)
        world.bind_to(model)
        world.add_offlattice(base, self.offlattice)

        world.add_molecules(species, 5)
        self.assertEqual(5, world.num_molecules(species))

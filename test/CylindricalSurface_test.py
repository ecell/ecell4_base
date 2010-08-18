#!/usr/bin/env python

import unittest
from egfrd import *
import _gfrd
import myrandom

class CylindricalSurfaceTestCase(unittest.TestCase):
    def setUp(self):
        self.radius = 1
        self.L = 10
        self.world = _gfrd.World(self.L, 3)

        # Random values.
        self.r1 = 1.6
        self.r2 = 2.1
        self.r3 = 2.7
        self.r4 = 4.2

    def cylinder_at_position(self, z):
        return create_cylindrical_surface('d',
                                          [5, 5, z],
                                          self.radius,
                                          [0, 0, 1],
                                          self.L)

    def test_random_positions(self):
        # For a cylinder along z-axis, position z should not matter.
        for z in range(0, self.L):
            d = self.cylinder_at_position(z)

            positions = []
            for i in range(100):
                position = d.random_position(myrandom.rng)
                position = apply_boundary(position, self.world.world_size)
                positions.append(position)

            average_position = numpy.average(positions, 0)
            assert average_position[0] == 5
            assert average_position[1] == 5
            assert 4 < average_position[2] < 6

    def test_projected_point(self):
        d = self.cylinder_at_position(self.r1)
        assert (d.projected_point([2, 2, 2]) == [5, 5, 2]).all()

    def test_distance_to_cylinder(self):
        d = self.cylinder_at_position(self.r1)

        assert self.world.distance(d.shape, [7, 5, self.r2]) == 1

        # Inside.
        assert self.world.distance(d.shape, [5, 5, self.r3]) == -1

if __name__ == "__main__":
    unittest.main()


#!/usr/bin/env python

import unittest
from egfrd import *
import _gfrd
import myrandom

class PlanarSurfaceTestCase(unittest.TestCase):
    def setUp(self):
        self.radius = 1
        self.L = 10
        self.world = _gfrd.World(self.L, 3)

        # Random values.
        self.r1 = 1.6
        self.r2 = 2.1
        self.r3 = 2.7
        self.r4 = 4.2

    def membrane_at_position(self, x, y):
        return create_planar_surface('m1',
                                     [x, y, 5],
                                     [1, 0, 0],
                                     [0, 1, 0],
                                     self.L,
                                     self.L)

    def test_random_positions(self):
        # For a plane in x-y, positions x and y should not matter.
        for x in range(0, self.L):
            for y in range(0, self.L):
                m1 = self.membrane_at_position(x, y)

                positions = []
                for i in range(100):
                    position = m1.random_position(myrandom.rng)
                    position = apply_boundary(position, self.world.world_size)
                    positions.append(position)

                average_position = numpy.average(positions, 0)
                assert 4 < average_position[0] < 6
                assert 4 < average_position[1] < 6
                assert average_position[2] == 5

    def test_projected_point(self):
        m = self.membrane_at_position(self.r1, self.r2)
        assert (m.projected_point([2, 2, 2]) == [2, 2, 5]).all()

    def test_distance_to_plane(self):
        m1 = self.membrane_at_position(self.r1, self.r2)

        assert self.world.distance(m1.shape, [self.r3, self.r4, 2]) == -3
        assert self.world.distance(m1.shape, [self.r3, self.r4, 8]) == 3

if __name__ == "__main__":
    unittest.main()

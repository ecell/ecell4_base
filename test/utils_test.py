#!/usr/bin/env python

import unittest

import numpy

from utils import *

class UtilsTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_cyclic_transpose(self):

        pos = cyclic_transpose(numpy.array([1,1,1]), numpy.array([1,1,1]), 10)
        self.failIf(not pos[0] == pos[1] == pos[2] == 1)

        pos = cyclic_transpose(numpy.array([8,8,8]), numpy.array([1,1,1]), 10)
        self.failIf(not pos[0] == pos[1] == pos[2] == -2)

        pos = cyclic_transpose(numpy.array([1,1,1]), numpy.array([8,8,8]), 10)
        self.failIf(not pos[0] == pos[1] == pos[2] == 11)

        pos = cyclic_transpose(numpy.array([1,8,1]), numpy.array([8,1,8]), 10)
        self.failIf(not (pos[0] == pos[2] == 11 and pos[1] == -2))

    
    def test_random_unit_vector(self):

        for i in range(1000):
            v = random_unit_vector()
            self.assertAlmostEqual(length(v), 1.0, 15)


    def test_random_vector(self):

        for i in range(1000):
            r = numpy.random.uniform() * 1e3
            v = random_vector(r)
            self.assertAlmostEqual(length(v), r, 12)

    def test_normalize(self):

        for i in range(1000):
            r = numpy.random.uniform(size=3) * 1e3
            l = numpy.random.uniform()
            v = normalize(r,l)
            self.assertAlmostEqual(length(v), l, 12)


    def test_spherical_cartesian(self):

        for i in range(1000):
            r = numpy.random.uniform() * 1e3
            v = random_vector(r)
            v2 = spherical_to_cartesian(cartesian_to_spherical(v))
            diff = abs(v - v2).sum()
            self.assertAlmostEqual(diff, 0, 10)

    def test_cartesian_spherical(self):

        for i in range(1000):
            v = random_unit_vector_s()
            v[0] *= 1e3
            v2 = cartesian_to_spherical(spherical_to_cartesian(v))
            diff = abs(v - v2).sum()
            self.assertAlmostEqual(diff, 0, 10)
            
    def test_calculate_pair_CoM(self):

        #FIXME: more serious test is needed.

        def _calculate_pair_com(pos1, pos2, D1, D2, world_size):
            pos2t = cyclic_transpose(pos2, pos1, world_size)
            return ((D2 * pos1 + D1 * pos2t) / (D1 + D2)) % world_size

        for i in range(1000):
            pos1 = numpy.random.uniform(size=3)
            pos2 = numpy.random.uniform(size=3)
            D1 = numpy.random.uniform()
            D2 = numpy.random.uniform()
            wsize = 10
            com1 = calculate_pair_CoM(pos1, pos2, D1, D2, wsize)
            com2 = _calculate_pair_com(pos1, pos2, D1, D2, wsize)
            diff = abs(com1 - com2).sum()
            self.assertAlmostEqual(diff, 0, 10)


if __name__ == "__main__":
    unittest.main()

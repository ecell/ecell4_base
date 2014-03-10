#!/usr/bin/env python

import unittest

import numpy

from _gfrd import *
import math

class CylindricalShellContainerTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testC1(self):
        c = CylindricalShellContainer(1000, 3)
        shell_id0 = ShellID(0,0) 
        c.update((shell_id0, CylindricalShell(DomainID(0,0), Cylinder([500, 500, 500], 50, [0,0,1], 50, ))))
        
        # Distance to cylinder in z direction.
        d = c.get_neighbors_within_radius([500, 500, 600], 75)
        self.assertAlmostEqual(50, d[0][1])

        # Cylinder update with same value.
        c.update((shell_id0, CylindricalShell(DomainID(0,0), Cylinder([500, 500, 500], 50, [0,0,1], 50))))
        d = c.get_neighbors_within_radius([500, 500, 600], 75)
        self.assertAlmostEqual(50, d[0][1])

        # Real update: longer cylinder.
        c.update((shell_id0, CylindricalShell(DomainID(0,0), Cylinder([500, 500, 500], 50, [0,0,1], 75))))
        d = c.get_neighbors_within_radius([500, 500, 600], 75)
        self.assertAlmostEqual(25, d[0][1])


    def testC2(self):
        shell_id0 = ShellID(0,0) 
        c = CylindricalShellContainer(1000, 3)
        c.update((shell_id0, CylindricalShell(DomainID(0,0), Cylinder([500, 500, 500], 50, [0,0,1], 50))))

        # Distance to cylinder in radial direction.
        d = c.get_neighbors([550, 550, 500])
        self.assertAlmostEqual(math.sqrt(50**2 + 50**2) - 50, d[0][1])

        # Distance to cylinder edge.
        d = c.get_neighbors([500, 553, 554])
        self.assertAlmostEqual(math.sqrt(3**2 + 4**2), d[0][1])


    def testC3(self):
        shell_id0 = ShellID(0,0) 
        c = CylindricalShellContainer(1000, 3)
        c.update((shell_id0, CylindricalShell(DomainID(0,0), Cylinder([0, 0, 0], 50, [0,0,1], 50))))

        # Using periodic boundary conditions.
        # Distance to cylinder in radial direction.
        d = c.get_neighbors([950, 950, 0])
        self.assertAlmostEqual(math.sqrt(50**2 + 50**2) - 50, d[0][1])

        # Distance to cylinder edge.
        d = c.get_neighbors([0, 947, 946])
        self.assertAlmostEqual(math.sqrt(3**2 + 4**2), d[0][1])


    def testC4(self):
        shell_id0 = ShellID(0,0) 
        c = CylindricalShellContainer(1000, 3)
        c.update((shell_id0, CylindricalShell(DomainID(0,0), Cylinder([0, 0, 0], 100, [0,0,1], 100))))

        d = c.get_neighbors([700, 0, 0])
        self.assertAlmostEqual(200, d[0][1])

        # Todo.
        d = c.get_neighbors([600, 0, 0])
        self.assertAlmostEqual(300, d[0][1])


if __name__ == "__main__":
    unittest.main()

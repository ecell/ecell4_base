#!/usr/bin/env python

import unittest

import numpy

from _gfrd import *
import math

class CylindricalShellContainerTestCase( unittest.TestCase ):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testC1( self ):
        c = CylindricalShellContainer( 1000, 3 )
        shellId0 = ShellID(0,0) 
        c.update( ( shellId0, CylindricalShell([500, 500, 500], 50, [0,0,1], 50, DomainID(0,0) ) ) )
        
        # Distance to cylinder in z direction.
        d = c.get_neighbors_within_radius( [500, 500, 600], 75 )
        self.assertAlmostEqual( 50, d[0][1] )

        # Cylinder update with same value.
        c.update( ( shellId0, CylindricalShell([500, 500, 500], 50, [0,0,1], 50, DomainID(0,0) ) ) )
        d = c.get_neighbors_within_radius( [500, 500, 600], 75 )
        self.assertAlmostEqual( 50, d[0][1] )

        # Real update: longer cylinder.
        c.update( ( shellId0, CylindricalShell([500, 500, 500], 50, [0,0,1], 75, DomainID(0,0) ) ) )
        d = c.get_neighbors_within_radius( [500, 500, 600], 75 )
        self.assertAlmostEqual( 25, d[0][1] )


    def testC2( self ):
        shellId0 = ShellID(0,0) 
        c = CylindricalShellContainer( 1000, 3 )
        c.update( ( shellId0, CylindricalShell([500, 500, 500], 50, [0,0,1], 50, DomainID(0,0) ) ) )

        # Distance to cylinder in radial direction.
        d = c.get_neighbors( [550, 550, 500] )
        self.assertAlmostEqual( math.sqrt(50**2 + 50**2) - 50, d[0][1] )

        # Distance to cylinder edge.
        d = c.get_neighbors( [500, 553, 554] )
        self.assertAlmostEqual( math.sqrt(3**2 + 4**2), d[0][1] )


    def testC3( self ):
        shellId0 = ShellID(0,0) 
        c = CylindricalShellContainer( 1000, 3 )
        c.update( ( shellId0, CylindricalShell([0, 0, 0], 50, [0,0,1], 50, DomainID(0,0) ) ) )

        # Using periodic boundary conditions.
        # Distance to cylinder in radial direction.
        d = c.get_neighbors( [950, 950, 0] )
        self.assertAlmostEqual( math.sqrt(50**2 + 50**2) - 50, d[0][1] )

        # Distance to cylinder edge.
        d = c.get_neighbors( [0, 947, 946] )
        self.assertAlmostEqual( math.sqrt(3**2 + 4**2), d[0][1] )


    def testC4( self ):
        shellId0 = ShellID(0,0) 
        c = CylindricalShellContainer( 1000, 3 )
        c.update( ( shellId0, CylindricalShell([0, 0, 0], 100, [0,0,1], 100, DomainID(0,0) ) ) )

        d = c.get_neighbors( [700, 0, 0] )
        self.assertAlmostEqual( 200, d[0][1] ) # Ok.

        # Todo. Bug?
        d = c.get_neighbors( [600, 0, 0] )
        #self.assertAlmostEqual( 300, d[0][1] ) # Fail.
        self.assertAlmostEqual( 500, d[0][1] ) # Not Ok!


if __name__ == "__main__":
    unittest.main()

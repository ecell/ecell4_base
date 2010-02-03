#!/usr/bin/env python

import unittest

import numpy

from _gfrd import *
import math

class SphericalShellContainerTestCase( unittest.TestCase ):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testS1( self ):
        c = SphericalShellContainer(1.0, 10)
        shellId0 = ShellID(0, 0)
        shellId1 = ShellID(0, 1)

        self.assertEqual( 10, c.matrix_size )
        self.assertAlmostEqual( 0.1, c.cell_size )
        self.assertEqual( 0, len(c) )

        c.update( ( shellId0, SphericalShell( [0.5, 0.3, 0.2], 0.1, DomainID(0, 0) ) ) )
        self.assertEqual( 1, len(c) )

        c.update( ( shellId1, SphericalShell( [0.0, 0.3, 0.9], 0.1, DomainID(0, 1) ) ) )
        self.assertEqual( 2, len(c) )

        self.assertAlmostEqual( c[shellId1].position[0], 0.0 )
        self.assertAlmostEqual( c[shellId1].position[1], 0.3 )
        self.assertAlmostEqual( c[shellId1].position[2], 0.9 )

        a = c.get_neighbors_within_radius( [0.45, 0.23, 0.13], 0.02 )
        # Distance to shell 0 is about 0.01 (should be found).
        # Distance to shell 1 is about 0.41 (should not be found).
        self.assertEqual( 1, len(a) )


    def testS2( self ):
        c = SphericalShellContainer( 1000, 3 )
        shellId0 = ShellID(0,0) 
        c.update( ( shellId0, SphericalShell( [500, 500, 500], 50, DomainID(0, 0) ) ) )
        
        # Find neighbors.
        d = c.get_neighbors_within_radius( [500, 500, 600], 75 )
        self.assertAlmostEqual( 50, d[0][1] )

        # Update with same value.
        # Returns false, but works fine.
        c.update( ( shellId0, SphericalShell( [500, 500, 500], 50, DomainID(0, 0) ) ) )
        d = c.get_neighbors_within_radius( [500, 500, 600], 75 )
        self.assertAlmostEqual( 50, d[0][1] )

        # Now a real update.
        # Returns false, but works fine.
        c.update( ( shellId0, SphericalShell( [500, 500, 500], 75, DomainID(0, 0) ) ) )
        d = c.get_neighbors_within_radius( [500, 500, 600], 100 )
        self.assertAlmostEqual( 25, d[0][1] )


    def testS3( self ):
        c = SphericalShellContainer( 1000, 3 )

        # A sphere at x=0 with radius=300.
        shellId0 = ShellID(0,0) 
        c.update( ( shellId0, SphericalShell( [0, 500, 500], 300, DomainID(0, 0) ) ) )

        # Distance to sphere from x=670 should be 30.
        d = c.get_neighbors_cyclic( [670, 500, 500] )
        self.assertAlmostEqual( 30, d[0][1] ) # Ok.

        # Distance to sphere from x=660 should be 40.
        d = c.get_neighbors_cyclic( [660, 500, 500] )
        self.assertAlmostEqual( 40, d[0][1] )


if __name__ == "__main__":
    unittest.main()

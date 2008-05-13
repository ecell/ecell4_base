#!/usr/bin/env python

import unittest

from numpy import ndarray

from object_matrix import *

class object_matrixTestCase( unittest.TestCase ):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test1( self ):

        c = ObjectContainer(1.0, 10)
        self.assertEqual( 10, c.matrix_size )
        self.assertEqual( 0.1, c.cell_size )
        self.assertEqual( 0, len(c) )
        c[0] = Sphere((0.5, 0.3, 0.2), 0.1)
        self.assertEqual( 1, len(c) )
        self.assertEqual( None, c[1] )
        tmp = ndarray(shape = (3, ))
        tmp[0] = 0.0
        tmp[1] = 0.3
        tmp[2] = 0.9
        c[1] = Sphere(tmp, 0.1)
        self.assertEqual( 2, len(c) )
        self.assertAlmostEqual( c[1].x, 0.0 )
        self.assertAlmostEqual( c[1].y, 0.3 )
        self.assertAlmostEqual( c[1].z, 0.9 )

        self.assertEqual( None, c[2] )
        self.failIf( not isinstance(c[0], SphereRef) )
        a = c.neighbors_array(Sphere((0.45, 0.23, 0.13), 0.09))
        self.assertEqual( 1, len(a[0]) )
        self.assertEqual( 1, len(a[1]) )
        self.assertEqual( 0, a[0][0].id )
        self.assertAlmostEqual( a[0][0].x, c[0].x )
        self.assertAlmostEqual( a[0][0].y, c[0].y )
        self.assertAlmostEqual( a[0][0].z, c[0].z )




if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python

import unittest

import numpy

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
        self.assertEqual( 0, c.size() )
        c.insert( 0, numpy.array([0.5, 0.3, 0.2]), 0.1 )
        self.assertEqual( 1, c.size() )
        #self.assertEqual( None, c.get( 1 ) )
        c.insert( 1, numpy.array([0.0,0.3,0.9]), 0.1)
        self.assertEqual( 2, c.size() )

        self.assertAlmostEqual( c.get(1)[0][0], 0.0 )
        self.assertAlmostEqual( c.get(1)[0][1], 0.3 )
        self.assertAlmostEqual( c.get(1)[0][2], 0.9 )

        #self.assertEqual( None, c[2] )

        a = c.neighbors_array( numpy.array([0.45, 0.23, 0.13]), 0.09 )
        self.assertEqual( 1, len(a[0]) )
        self.assertEqual( 1, len(a[1]) )
        self.assertEqual( 0, a[0][0] )





if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python

import unittest

import numpy

from utils import *

class UtilsTestCase( unittest.TestCase ):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_CyclicTranspose( self ):

        pos = cyclicTranspose( numpy.array([1,1,1]), numpy.array([1,1,1]), 10 )
        self.failIf( not pos[0] == pos[1] == pos[2] == 1 )

        pos = cyclicTranspose( numpy.array([8,8,8]), numpy.array([1,1,1]), 10 )
        self.failIf( not pos[0] == pos[1] == pos[2] == -2 )

        pos = cyclicTranspose( numpy.array([1,1,1]), numpy.array([8,8,8]), 10 )
        self.failIf( not pos[0] == pos[1] == pos[2] == 1 )

    
    def test_distanceSq_Cyclic( self ):

        distance = distanceSq_Cyclic( numpy.array([1.,1.,1.]),
                                      numpy.array([3.,3.,3.]), 10 )
        self.assertAlmostEqual( numpy.sqrt( 12. ), numpy.sqrt( distance ) )

        distance = distanceSq_Cyclic( numpy.array([1.,1.,1.]),
                                      numpy.array([9.,9.,9.]), 10 )
        self.assertAlmostEqual( numpy.sqrt( 12. ), numpy.sqrt( distance ) )

        distance = distanceSq_Cyclic( numpy.array([9.,9.,9.]),
                                      numpy.array([1.,1.,1.]), 10 )
        self.assertAlmostEqual( numpy.sqrt( 12. ), numpy.sqrt( distance ) )

        distance = distanceSq_Cyclic( numpy.array([9.,1.,9.]),
                                      numpy.array([1.,9.,1.]), 10 )
        self.assertAlmostEqual( numpy.sqrt( 12. ), numpy.sqrt( distance ) )

    def test_distanceSqArray_Cyclic( self ):

        distance = distanceSqArray_Cyclic( numpy.array([1.,1.,1.]),
                                           ( numpy.array([3.,3.,3.]),
                                             numpy.array([9.,9.,9.]) ), 10 )

        self.assertAlmostEqual( numpy.sqrt( 12. ), numpy.sqrt( distance[0] ) )
        self.assertAlmostEqual( numpy.sqrt( 12. ), numpy.sqrt( distance[1] ) )

        distance = distanceSqArray_Cyclic( numpy.array([9.,9.,9.]),
                                           ( numpy.array([1.,1.,1.]),
                                             numpy.array([7.,7.,7.]) ), 10 )
        self.assertAlmostEqual( numpy.sqrt( 12. ), numpy.sqrt( distance[0] ) )
        self.assertAlmostEqual( numpy.sqrt( 12. ), numpy.sqrt( distance[1] ) )

        distance = distanceSqArray_Cyclic( numpy.array([9.,1.,9.]),
                                           ( numpy.array([1.,9.,1.]),
                                             numpy.array([7.,3.,7.]) ), 10 )
        self.assertAlmostEqual( numpy.sqrt( 12. ), numpy.sqrt( distance[0] ) )
        self.assertAlmostEqual( numpy.sqrt( 12. ), numpy.sqrt( distance[1] ) )




if __name__ == "__main__":
    unittest.main()

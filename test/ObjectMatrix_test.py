#!/usr/bin/env python

import unittest

import numpy

import ObjectMatrix as mod


class Obj( object ):
    def __init__( self, pos, radius ):
        self.pos = numpy.array( pos )
        self.radius = radius


class ObjectMatrixTestCase( unittest.TestCase ):

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def testInstantiation(self):
        m = mod.ObjectMatrix()
        self.failIf( m is None )

    def testEmptyState(self):
        m = mod.ObjectMatrix()
        self.failIf( m.size != 0 )

    def testSingle(self):
        m = mod.ObjectMatrix()
        m.setWorldSize( 1.0 )

        o1 = Obj( [0,0,.1], .1 )
        m.add( o1 )
        self.assertEqual( m.size, 1 )

        n, d = m.getNeighbors( numpy.array( [0,0,0] ) )

        self.failIf( len( n ) != 1 and len( d ) != 1 )
        self.failIf( n[0] != o1 )
        self.failIf( d[0] != 0 )

    def testTwo(self):
        m = mod.ObjectMatrix()

        o1 = Obj( [0,0,.1], .1 )
        o2 = Obj( [0,.3,0], .1 )
        m.add( o1 )
        m.add( o2 )
        self.assertEqual( m.size, 2 )

        n, d = m.getNeighbors( numpy.array( [0,0,0] ) )

        self.failIf( len( n ) != 2 and len( d ) != 2 )
        self.failIf( n[0] != o1 and n[1] != o2 )
        self.failIf( d[0] != 0 and d[1] != .2 )



if __name__ == "__main__":
    unittest.main()

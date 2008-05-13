#!/usr/bin/env python

import unittest

import numpy

import cObjectMatrix as mod


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
        m.add( o1, o1.pos, o1.radius )
        self.assertEqual( m.size, 1 )

        n, d = m.getNeighborsCyclic( numpy.array( [0,0,0] ) )
        print d
        self.failIf( len( n ) != 1 or len( d ) != 1 )
        self.failIf( n[0] != o1 )
        self.failIf( d[0] != 0 )

    def testTwo(self):
        m = mod.ObjectMatrix()
        m.setWorldSize( 1.0 )

        o1 = Obj( [0,0,.1], .1 )
        o2 = Obj( [0,.3,0], .1 )
        m.add( o1, o1.pos, o1.radius )
        m.add( o2, o2.pos, o2.radius )
        self.assertEqual( m.size, 2 )

        n, d = m.getNeighborsCyclic( numpy.array( [0,0,0] ) )

        self.failIf( len( n ) != 2 or len( d ) != 2 )
        self.failIf( n[0] != o1 or n[1] != o2 )
        self.assertAlmostEqual( 0, d[0] )
        self.assertAlmostEqual( .2, d[1] )


    def testCyclic(self):
        m = mod.ObjectMatrix()
        m.setWorldSize( 1.0 )

        o1 = Obj( [0,0,.1], .1 )
        o2 = Obj( [0,.8,0], .1 )
        m.add( o1, o1.pos, o1.radius )
        m.add( o2, o2.pos, o2.radius )
        self.assertEqual( m.size, 2 )

        n, d = m.getNeighborsCyclic( numpy.array( [0,0,0] ) )

        self.failIf( len( n ) != 2 or len( d ) != 2 )
        self.assertAlmostEqual( 0, d[0] )
        self.assertEqual( o1, n[0] )
        self.assertAlmostEqual( .1, d[1] )
        self.assertEqual( o2, n[1] )


    def testRemove(self):
        m = mod.ObjectMatrix()
        m.setWorldSize( 1.0 )

        o1 = Obj( [0,0,.1], .1 )
        o2 = Obj( [0,.8,0], .1 )
        o3 = Obj( [0,.9,.9], .1 )
        m.add( o1, o1.pos, o1.radius )
        self.assertEqual( m.size, 1 )
        m.remove( o1 )
        self.assertEqual( m.size, 0 )

        m.add( o1, o1.pos, o1.radius )
        m.add( o2, o2.pos, o2.radius )
        m.add( o3, o3.pos, o3.radius )
        self.assertEqual( m.size, 3 )

        m.remove( o2 )
        self.assertEqual( m.size, 2 )
        n, d = m.getNeighborsCyclic( numpy.array( [0,0,0] ) )
        self.failIf( len( n ) != 2 or n[0] != o1 )

        m.remove( o1 )
        self.assertEqual( m.size, 1 )

        m.remove( o3 )
        self.assertEqual( m.size, 0 )


    def testUpdate(self):
        m = mod.ObjectMatrix()
        m.setWorldSize( 1.0 )

        o1 = Obj( [0,0,.1], .1 )
        o2 = Obj( [0,.8,0], .1 )
        m.add( o1, o1.pos, o1.radius )
        m.add( o2, o2.pos, o2.radius )
        self.assertEqual( m.size, 2 )


        m.update( o1, [.5,.5,.5],.2 )
        n, d = m.getNeighborsCyclic( numpy.array( [.5,.5,.7] ) )
        self.failIf( len( n ) != 2 or n[0] != o1 )
        self.assertAlmostEqual( 0, d[0] )



if __name__ == "__main__":
    unittest.main()

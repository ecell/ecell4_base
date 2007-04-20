#!/usr/bin/env python

import unittest

import _gfrd as mod

class FirstPassageGreensFunctionTestCase( unittest.TestCase ):

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_instantiation( self ):
        D = 1e-12
        gf = mod.FirstPassageGreensFunction( D )
        self.failIf( gf == None )

    def test_drawR_smallt( self ):
        D = 1e-12
        a = 1e-8
        gf = mod.FirstPassageGreensFunction( D )
        t = gf.drawTime( 0.5, a )
        t *= 1e-4
        r = gf.drawR( 0.5, t, a )
        self.failIf( r <= 0.0 )
        self.failIf( r > a )

    def test_drawR_large_t( self ):
        D = 1e-12
        a = 1e-6
        gf = mod.FirstPassageGreensFunction( D )
        t = gf.drawTime( 1.0 - 1e-16, a )
        print 't', t
        r = gf.drawR( 0.5, t, a )
        print r
        self.failIf( r <= 0.0 )
        self.failIf( r > a )


if __name__ == "__main__":
    unittest.main()

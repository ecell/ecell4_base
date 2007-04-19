#!/usr/bin/env python

import unittest

import _gfrd as mod

import numpy


class FirstPassagePairGreensFunctionTestCase( unittest.TestCase ):

    def setUp( self ):
        self.N8 = 0.99999999

    def tearDown( self ):
        pass
    
    def testInstantiation( self ):
        D = 1e-12
        kf = 1e8
        Sigma = 1e-8
        a = 1e-7

        gf = mod.FirstPassagePairGreensFunction( D, kf, Sigma )
        self.failIf( gf == None )
        gf.seta( a )


    def testDrawTime( self ):
        D = 1e-12
        kf = 1e8
        Sigma = 1e-8
        a = 1e-7
        r0 = 5e-8
        
        gf = mod.FirstPassagePairGreensFunction( D, kf, Sigma )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        self.failIf( t <= 0.0 or t >= numpy.inf )

        t = gf.drawTime( 0.0, r0 )
        self.failIf( t < 0.0 or t >= numpy.inf )

        t = gf.drawTime( 1.0, r0 )
        self.failIf( t <= 0.0 or t >= numpy.inf )

    def testDrawEventType( self ):
        D = 1e-12
        kf = 1e8
        Sigma = 1e-8
        a = 1e-7
        r0 = 5e-8
        
        gf = mod.FirstPassagePairGreensFunction( D, kf, Sigma )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        eventType = gf.drawEventType( 0.5, r0, t )
        self.failIf( eventType != 0 and eventType != 1 and eventType != 2 )

        eventType = gf.drawEventType( 0.0, r0, t )
        self.failIf( eventType != 0 and eventType != 1 and eventType != 2 )

        eventType = gf.drawEventType( 1.0, r0, t )
        self.failIf( eventType != 0 and eventType != 1 and eventType != 2 )


    def testDrawR( self ):
        D = 1e-12
        kf = 1e8
        Sigma = 1e-8
        a = 1e-7
        r0 = 5e-8
        
        gf = mod.FirstPassagePairGreensFunction( D, kf, Sigma )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        eventType = gf.drawEventType( 0.5, r0, t )
        self.failIf( eventType != 0 and eventType != 1 and eventType != 2 )

        r = gf.drawR( 0.5, r0, t )
        self.failIf( r < Sigma or r > a )

        r = gf.drawR( 0.0, r0, t )
        self.failIf( r < Sigma or r > a )
        self.assertEqual( r, Sigma )

        r = gf.drawR( self.N8, r0, t )
        self.failIf( r < Sigma or r > a )
        #self.failIf( r != a )


    def testDrawTheta( self ):
        D = 1e-12
        kf = 1e8
        Sigma = 1e-8
        a = 1e-7
        r0 = 5e-8
        
        gf = mod.FirstPassagePairGreensFunction( D, kf, Sigma )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        eventType = gf.drawEventType( 0.5, r0, t )
        r = gf.drawR( 0.5, r0, t )

        theta = gf.drawTheta( 0.5, r, r0, t )
        self.failIf( theta < 0.0 or theta > 2 * numpy.pi )

        theta = gf.drawTheta( 0.0, r, r0, t )
        self.failIf( theta < 0.0 or theta > 2 * numpy.pi )

        theta = gf.drawTheta( 1.0, r, r0, t )
        self.failIf( theta < 0.0 or theta > 2 * numpy.pi )



    def testAlpha0( self ):

        D = 1e-12
        Sigma = 1e-8
        kf = 1e-18
        
        a = 2e-7
        
        gf = mod.FirstPassagePairGreensFunction( D, kf, Sigma )
        gf.seta( a )
        maxerror = 0
        
        for i in range(1000):
            alpha = gf.alpha0_i( i )
            error = abs( gf.f_alpha0( alpha ) )
            maxerror = max( error, maxerror )

        self.failIf( abs( maxerror ) > 1e-8 )


        
if __name__ == "__main__":
    unittest.main()

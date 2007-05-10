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
        kf = 1e-8
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

    def testDrawTime_a_equal_sigma( self ):
        D = 1e-12
        kf = 1e-8
        Sigma = 1e-8
        a = Sigma
        r0 = a
        
        gf = mod.FirstPassagePairGreensFunction( D, kf, Sigma )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        self.assertEqual( 0.0, t )

    def testDrawTime_r0_equal_a( self ):
        D = 1e-12
        kf = 1e-8
        Sigma = 1e-8
        a = 1e-7
        r0 = a
        
        gf = mod.FirstPassagePairGreensFunction( D, kf, Sigma )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        self.assertEqual( 0.0, t )


    def testDrawEventType( self ):
        D = 1e-12
        kf = 1e-8
        Sigma = 1e-8
        a = 1e-7
        r0 = 5e-8

        gf = mod.FirstPassagePairGreensFunction( D, kf, Sigma )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        eventType = gf.drawEventType( 0.5, r0, t )
        self.failIf( eventType != 0 and eventType != 1 and eventType != 2 )

        eventType = gf.drawEventType( 0.0, r0, t )
        self.assertEqual( eventType, 0 )

        eventType = gf.drawEventType( 1.0, r0, t )
        self.assertEqual( eventType, 1 )


    def testDrawR( self ):
        D = 1e-12
        kf = 1e-8
        Sigma = 1e-8
        a = 1e-7
        r0 = 2e-8
        
        gf = mod.FirstPassagePairGreensFunction( D, kf, Sigma )
        gf.seta( a )

        t = 1e-3

        r = gf.drawR( 0.5, r0, t )
        self.failIf( r < Sigma or r > a )

        r1 = gf.drawR( 0.0, r0, t )
        r2 = gf.drawR( 1.0, r0, t )

        self.failIf( r1 < Sigma or r1 > a )
        self.failIf( r2 < Sigma or r2 > a )

        self.failIf( abs( r1 - Sigma ) > 1e-15 )
        self.failIf( abs( r2 - a ) > 1e-15 )


    def testDrawR_zerot( self ):
        D = 1e-12
        kf = 1e-8
        Sigma = 1e-8
        a = 1e-7
        r0 = 2e-8
        
        gf = mod.FirstPassagePairGreensFunction( D, kf, Sigma )
        gf.seta( a )

        t = 0.0

        r = gf.drawR( 0.5, r0, t )
        self.assertEqual( r0, r )

    def testDrawTheta( self ):
        D = 1e-12
        kf = 1e-8
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


    def testDrawTheta_zerot( self ):
        D = 1e-12
        kf = 1e-8
        Sigma = 1e-8
        a = 1e-7
        r = 5e-8
        r0 = 5e-8
        
        gf = mod.FirstPassagePairGreensFunction( D, kf, Sigma )
        gf.seta( a )

        t = 0.0
        theta = gf.drawTheta( 0.5, r0, r0, t )
        self.assertEqual( 0.0, theta )



    def testAlpha0( self ):

        D = 1e-12
        Sigma = 1e-8
        kf = 1e-8
        
        a = 1e-7
        
        gf = mod.FirstPassagePairGreensFunction( D, kf, Sigma )
        gf.seta( a )
        maxerror = 0.0
        
        for i in range(100):
            alpha = gf.alpha0_i( i )
            error = abs( gf.f_alpha0( alpha ) / alpha )
            #print error/alpha, gf.f_alpha0( alpha*1.1 )/alpha
            maxerror = max( error, maxerror )

        self.failIf( abs( maxerror ) > 1e-10 )

'''
    def testAlphan( self ):

        D = 1e-12
        Sigma = 1e-8
        kf = 1e-18
        
        a = 2e-7
        
        gf = mod.FirstPassagePairGreensFunction( D, kf, Sigma )
        gf.seta( a )
        maxerror = 0
        
        for n in range(100):
            for i in range(1000):
                alpha = gf.alpha_i( n, i )
                error = abs( gf.f_alpha0( alpha ) )
                maxerror = max( error, maxerror )

        self.failIf( abs( maxerror ) > 1e-8 )
'''

        
if __name__ == "__main__":
    unittest.main()

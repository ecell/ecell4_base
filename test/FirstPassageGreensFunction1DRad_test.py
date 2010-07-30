#!/usr/bin/env python

__author__    = 'Laurens Bossen, Thomas Sokolowski'
__copyright__ = ''


import unittest

import _greens_functions as mod

import numpy


class FirstPassageGreensFunction1DRadTestCase( unittest.TestCase ):

    def setUp( self ):
        pass

    def tearDown( self ):
        pass

    def test_Instantiation( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 2e-7
	r0 = L/2

        gf = mod.FirstPassageGreensFunction1DRad( D, kf, v, r0, L )
        self.failIf( gf == None )


    def test_DrawTime( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 2e-7
        r0 = 5e-8

        gf = mod.FirstPassageGreensFunction1DRad( D, kf, v, r0, L )

        t = gf.drawTime( 0.5 )
        self.failIf( t <= 0.0 or t >= numpy.inf )
        print " "
        print "FirstPassageGreensFunction1DRad_test.py : test_DrawTime : t =",t

        t = gf.drawTime( 0.0 )
        self.failIf( t < 0.0 or t >= numpy.inf )
        print "FirstPassageGreensFunction1DRad_test.py : test_DrawTime : t =",t

        t = gf.drawTime( 1 - 1e-16 )
        self.failIf( t <= 0.0 or t >= numpy.inf )
        print "FirstPassageGreensFunction1DRad_test.py : test_DrawTime : t =",t


    def test_DrawTime_a_equal_sigma( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 0
        r0 = L

        gf = mod.FirstPassageGreensFunction1DRad( D, kf, v, r0, L )
        t = gf.drawTime( 0.5 )
        self.assertEqual( 0.0, t )
        print " "
        print "FirstPassageGreensFunction1DRad_test.py : test_DrawTime_a_equal_sigma : t =",t


    def test_DrawTime_a_near_sigma( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 2e-14
        r0 = L/2

        gf = mod.FirstPassageGreensFunction1DRad( D, kf, v, r0, L )

        print "drawTime...a near sigma"
        t = gf.drawTime( 0.5 )
        print "done"
        self.failIf( t <= 0.0 or t >= numpy.inf )
        print " "
        print "FirstPassageGreensFunction1DRad_test.py : test_DrawTime_a_near_sigma : t =",t


    def test_DrawTime_r0_equal_a( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 2e-7
        r0 = L

        gf = mod.FirstPassageGreensFunction1DRad( D, kf, v, r0, L )

        t = gf.drawTime( 0.5 )
        self.assertEqual( 0.0, t )
        print " "
        print "FirstPassageGreensFunction1DRad_test.py : test_DrawTime_r0_equal_a : t =",t


    def test_DrawTime_r0_equal_sigma_kf_zero( self ):
        D = 1e-12
	v = -3e-8
        kf = 0.0 # note this
        L = 1e-7
        r0 = 0

        gf = mod.FirstPassageGreensFunction1DRad( D, kf, v, r0, L )

        t = gf.drawTime( 0.5 )				# This does not converge for negative v ~< -1e-16
        self.failIf( t < 0.0 or t >= numpy.inf )
        print " "
        print "FirstPassageGreensFunction1DRad_test.py : test_DrawTime_r0_equal_sigma_kf_zero : t =",t


    def no_test_DrawTime_r0_equal_sigma_kf_large( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 20e-7
        r0 = 1e-12

        gf = mod.FirstPassageGreensFunction1DRad( D, kf, v, r0, L )

        t = gf.drawTime( 0.5 )
        self.failIf( t < 0.0 or t >= numpy.inf )
        print " "
        print "FirstPassageGreensFunction1DRad_test.py : test_DrawTime_r0_equal_sigma_kf_large : t =",t


    def test_DrawEventType( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 2e-7
        r0 = L/2

        gf = mod.FirstPassageGreensFunction1DRad( D, kf, v, r0, L )

        t = gf.drawTime( 0.5 )
        eventType = gf.drawEventType( 0.5, t )
        self.failIf( eventType != 12 and eventType != 13 and eventType != 14 )
        print " "
        print "FirstPassageGreensFunction1DRad_test.py : test_DrawEventType : eventType =",eventType

        eventType = gf.drawEventType( 0.999999, t )
        self.assertEqual( eventType, 13 ) # ESCAPE

        eventType = gf.drawEventType( 0.0, t )
        self.assertEqual( eventType, 14 ) # REACTION


    def no_test_DrawEventType_smallt( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 2e-6
        r0 = L/2

        gf = mod.FirstPassageGreensFunction1DRad( D, kf, v, r0, L )

        t = gf.drawTime( 0.001 )

        eventType = gf.drawEventType( 0.5, t )
        self.failIf( eventType != 12 and eventType != 13 and eventType != 14 )
        print " "
        print "FirstPassageGreensFunction1DRad_test.py : test_DrawEventType_smallt : eventType =",eventType

        eventType = gf.drawEventType( 0.9999, t )
        self.assertEqual( eventType, 13 ) # ESCAPE

        eventType = gf.drawEventType( 0.0, t )
        self.assertEqual( eventType, 14 ) # REACTION


    def test_DrawR( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 2e-7
        r0 = L/2

        gf = mod.FirstPassageGreensFunction1DRad( D, kf, v, r0, L )

        t = 1e-3

        r = gf.drawR( 0.5, t )
        self.failIf( r < 0 or r > L )

        r1 = gf.drawR( 0.0, t )
        r2 = gf.drawR( 0.999999999999, t )

        self.failIf( r1 < 0 or r1 > L )
        self.failIf( r2 < 0 or r2 > L )

        self.assertAlmostEqual( r1, 0 )
        self.assertAlmostEqual( r2, L )
        print " "
        print "FirstPassageGreensFunction1DRad_test.py : test_DrawR : r1 =",r1
        print "FirstPassageGreensFunction1DRad_test.py : test_DrawR : r2 =",r2


    def test_DrawR_zerot( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 1e-7
        r0 = L/2

        gf = mod.FirstPassageGreensFunction1DRad( D, kf, v, r0, L )

        t = 0.0

        r = gf.drawR( 0.5, t )
        self.assertEqual( r0, r )
        print " "
        print "FirstPassageGreensFunction1DRad_test.py : test_DrawR_zerot : r =",r


    def test_DrawR_r0_equal_sigma( self ):
        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 2e-7
        r0 = 0

        t = 1e-3

        gf = mod.FirstPassageGreensFunction1DRad( D, kf, v, r0, L )

        r = gf.drawR( 0.5, t )
        self.failIf( r < 0 or r > L )


    def test_DrawR_squeezed( self ):

        D = 1e-12
	v = -3e-8
        kf = 1e-8
        L = 0.02e-8
	r0 = L/2

        gf = mod.FirstPassageGreensFunction1DRad( D, kf, v, r0, L )

        t = 1e-6
        r0 = 0
        gf.setr0 ( r0 )
        r = gf.drawR( 0.5, t )
        self.failIf( r < 0 or r > L )

        # near s
        r0 = 0.0001e-8
        gf.setr0 ( r0 )
        r = gf.drawR( 0.5, t )
        self.failIf( r < 0 or r > L )
        print " "
        print "FirstPassageGreensFunction1DRad_test.py : test_DrawR_squeezed : 1st drawn r =",r

        # near a
        r0 = L - 0.0001e-8
        gf.setr0 ( r0 )
        r = gf.drawR( 0.5, t )
        self.failIf( r < 0 or r > L )
        self.failIf( r < 0 or r > L )
        print "FirstPassageGreensFunction1DRad_test.py : test_DrawR_squeezed : 2nd drawn r =",r

if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python

__author__    = 'Koichi Takahashi <shafi@e-cell.org>'
__license__   = 'GPL'
__copyright__ = 'Copyright The Molecular Sciences Institute 2006-2008'


import unittest

import _gfrd as mod

import numpy


class PlainPairGreensFunctionTestCase( unittest.TestCase ):

    def setUp( self ):
        pass

    def tearDown( self ):
        pass
    
    def test_Instantiation( self ):
        D = 1e-12
        kf = 1e8
        sigma = 1e-8
        a = 1e-7

        gf = mod.PlainPairGreensFunction( D, kf, sigma )
        self.failIf( gf == None )


    def test_DrawTime( self ):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        r0 = 5e-8
        
        gf = mod.PlainPairGreensFunction( D, kf, sigma )

        t = gf.drawTime( 0.5, r0 )
        self.failIf( t <= 0.0 )

        t = gf.drawTime( 0.0, r0 )
        self.failIf( t < 0.0 )

        t = gf.drawTime( 1.0, r0 )
        self.failIf( t <= 0.0 )


    def test_DrawTime_r0_equal_sigma_kf_zero( self ):
        D = 1e-12
        kf = 0.0 # note this
        sigma = 1e-8
        r0 = sigma
        
        gf = mod.PlainPairGreensFunction( D, kf, sigma )

        t = gf.drawTime( 0.5, r0 )
        self.failIf( t < 0.0 or t >= numpy.inf )

    def test_DrawR( self ):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        r0 = 2e-8
        
        gf = mod.PlainPairGreensFunction( D, kf, sigma )

        t = 1e-3

        r = gf.drawR( 0.5, r0, t )
        self.failIf( r < sigma )

        r1 = gf.drawR( 0.0, r0, t )
        r2 = gf.drawR( 1.0, r0, t )

        self.failIf( r1 < sigma )
        self.failIf( r2 < sigma )

        self.failIf( abs( r1 - sigma ) > 1e-15 )


    def test_DrawR_zerot( self ):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        r0 = 2e-8
        
        gf = mod.PlainPairGreensFunction( D, kf, sigma )

        t = 0.0

        r = gf.drawR( 0.5, r0, t )
        self.assertEqual( r0, r )


    def test_DrawR_r0_equal_sigma( self ):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        r0 = sigma

        t = 1e-3
        
        gf = mod.PlainPairGreensFunction( D, kf, sigma )

        r = gf.drawR( 0.5, r0, t )
        self.failIf( r < sigma )


    def test_DrawTheta( self ):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        r0 = 5e-8
        
        gf = mod.PlainPairGreensFunction( D, kf, sigma )

        t = gf.drawTime( 0.5, r0 )
        r = gf.drawR( 0.5, r0, t )

        theta = gf.drawTheta( 0.5, r, r0, t )
        self.failIf( theta < 0.0 or theta > numpy.pi )

        theta = gf.drawTheta( 0.0, r, r0, t )
        self.failIf( theta < 0.0 or theta > numpy.pi )

        theta = gf.drawTheta( 1.0, r, r0, t )
        self.failIf( theta < 0.0 or theta > numpy.pi )


    def test_DrawTheta_zerot( self ):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        r = 5e-8
        r0 = 5e-8
        
        gf = mod.PlainPairGreensFunction( D, kf, sigma )

        t = 0.0
        theta = gf.drawTheta( 0.5, r0, r0, t )
        self.assertEqual( 0.0, theta )

    def test_DrawTheta_smallt( self ):

        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        r = 5e-8
        r0 = 5e-8
        
        gf = mod.PlainPairGreensFunction( D, kf, sigma )

        t = 1e-4  # well this is not *very* small..
        theta = gf.drawTheta( 0.5, r, r0, t )
        self.failIf( theta < 0.0 or theta > numpy.pi )


    def test_DrawTheta_r0_equal_sigma( self ):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        r0 = sigma

        t = 1e-3
        r = r0
        
        gf = mod.PlainPairGreensFunction( D, kf, sigma )

        theta = gf.drawTheta( 0.5, r, r0, t )
        self.failIf( theta < 0.0 or theta > numpy.pi )

    def test_ip_theta_is_int_p_theta( self ):

        import scipy.integrate

        D = 1e-12
        sigma = 1e-8
        kf = 1e-8

        t = 1e-2  #FIXME: smaller t should be fine
        r0 = 5e-8

        gf = mod.PlainPairGreensFunction( D, kf, sigma )
        r = r0

        ip = gf.ip_theta( 0.0, r, r0, t )
        self.assertEqual( 0.0, ip )
        
        resolution = 10
        for i in range( 1, resolution ):
            theta = i * numpy.pi / resolution 
            ip = gf.ip_theta( theta, r, r0, t )
            result = scipy.integrate.quad( gf.p_theta, 0.0, theta,
                                           args=( r, r0, t ) )
            np = result[0]
            #print theta, np, ip
            self.assertAlmostEqual( 0.0, (np-ip)/ip )


    def test_p_theta_never_negative( self ):

        D = 1e-12
        sigma = 1e-8
        kf = 1e-8

        # smaller t causes problem
        t = 1e-3
        r0 = 5e-8
        r = r0
        
        gf = mod.PlainPairGreensFunction( D, kf, sigma )

        pint = gf.ip_theta( numpy.pi, r, r0, t )

        pmin = 0.0
        resolution = 50
        for i in range( resolution ):
            theta = i * numpy.pi / resolution
            p = gf.p_theta( theta, r, r0, t ) / pint / resolution 
            pmin = min( pmin, p )
            #print 'theta: ', theta, '\tp: ', p
            
        self.failIf( pmin < 0.0, 'Negative p_theta; t= %g, %s'
                     % ( t, gf.dump() ) )


    def test_ip_theta_never_decrease( self ):

        D = 1e-12
        sigma = 1e-8
        kf = 1e-8

        # smaller t causes problem
        t = 1e-3
        r0 = 5e-8
        r = r0
        
        gf = mod.PlainPairGreensFunction( D, kf, sigma )

        pint_prev = 0.0

        resolution = 50
        for i in range( resolution ):
            theta = i * numpy.pi / resolution
            pint = gf.ip_theta( theta, r, r0, t )
            self.failIf( pint < pint_prev )
            pint_prev = pint



        
if __name__ == "__main__":
    unittest.main()

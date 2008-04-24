#!/usr/bin/env python

__author__    = 'Koichi Takahashi <shafi@e-cell.org>'
__license__   = 'GPL'
__copyright__ = 'Copyright The Molecular Sciences Institute 2006-2007'


import unittest

import _gfrd as mod

import math
import numpy


class FreeFunctionsTestCase( unittest.TestCase ):

    def setUp( self ):
        pass

    def tearDown( self ):
        pass

    def test_int_p_theta_free_is_ip_theta_free( self ):

        import scipy.integrate

        D = 1e-12
        t = 1e-5
        sigma = 1e-9
        r0 = 1e-9
        r = r0
        kf = 1e-18
        
        ip = mod.ip_theta_free( 0.0, r, r0, t, D )
        self.assertEqual( 0.0, ip )
        
        resolution = 10
        for i in range( 1, resolution ):
            theta = i * numpy.pi / resolution 
            ip = mod.ip_theta_free( theta, r, r0, t, D )
            result = scipy.integrate.quad( mod.p_theta_free, 0.0, theta,
                                           args=( r, r0, t, D ) )
            np = result[0]
            self.assertAlmostEqual( 0.0, (np-ip)/ip )


    def test_int_p_irr_is_S_irr( self ):

        import scipy.integrate

        D = 1e-12
        t = 1e-5
        sigma = 1e-9
        r0 = 1e-9
        kf = 1e-18
        

        for i in range( 1, 20 ):
            S = mod.S_irr( t, r0 * i, kf, D, sigma )
            result = scipy.integrate.quad( mod.p_irr, sigma, sigma * 1e3,
                                           args=( t, r0 * i, kf, D, sigma ) )
            ip = result[0]
            self.failIf( ip == 0 )
            self.assertAlmostEqual( 0.0, (S-ip)/ip )

    def test_int_g_bd_is_I_bd( self ):

        import scipy.integrate
        import math

        D = 1e-12
        t = 1e-6
        sigma = 1e-8
        r0 = 1e-9

        ibd = mod.I_bd( sigma, t, D )
        print ibd
        result = scipy.integrate.quad( mod.g_bd, sigma, 
                                       sigma + 6 * math.sqrt( 6 * D * t ),
                                       args=( sigma, t, D ) )
        igbd = result[0]
        print igbd
        self.failIf( ibd == 0 )
        self.assertAlmostEqual( 0.0, (ibd-igbd)/ibd )


    def test_int_g_bd_is_I_bd_smallt( self ):

        import scipy.integrate

        D = 1e-12
        t = 1e-20
        sigma = 1e-8
        r0 = 1e-9

        ibd = mod.I_bd( sigma, t, D )
        print ibd
        result = scipy.integrate.quad( mod.g_bd, sigma, sigma + 
                                       6 * math.sqrt( 6 * D * t ),
                                       args=( sigma, t, D ) )
        igbd = result[0]
        print igbd
        self.failIf( ibd == 0 )
        self.assertAlmostEqual( 0.0, (ibd-igbd)/ibd )


    def test_I_bd_r_large_is_I_bd( self ):

        D = 1e-12
        t = 1e-10
        sigma = 1e-8
        r0 = 1e-9

        ibd = mod.I_bd( sigma, t, D )
        ibdr = mod.I_bd_r( sigma + 6 * math.sqrt( 6 * D * t ), sigma, t, D )
        print ibd, ibdr

        self.assertAlmostEqual( 0.0, (ibd-ibdr)/ibd )


    def test_int_g_bd_is_I_bd_r( self ):

        import scipy.integrate
        import math

        D = 1e-12
        t = 1e-7
        sigma = 1e-8

        r_max = 6 * math.sqrt( 6 * D * t )

        ibd = mod.I_bd_r( sigma, sigma, t, D )
        self.failIf( ibd != 0.0 )

        N = 20
        for i in range( 1, N ):
            r = sigma + r_max / N * i
            ibd = mod.I_bd_r( r, sigma, t, D )
            result = scipy.integrate.quad( mod.g_bd, sigma, r,
                                           args=( sigma, t, D ) )
            igbd = result[0]

            self.failIf( ibd == 0 )
            self.assertAlmostEqual( 0.0, (ibd-igbd)/ibd )


    def test_drawR_gbd( self ):

        import scipy.integrate
        import math

        D = 1e-12
        t = 1e-10
        sigma = 1e-8

        r = mod.drawR_gbd( 0.0, sigma, t, D )
        self.assertEqual( r, sigma )

        r = mod.drawR_gbd( 0.5, sigma, t, D )
        self.failIf( r <= sigma )
        #print 'rr', r

        r = mod.drawR_gbd( 1.0, sigma, t, D )
        self.failIf( r <= sigma )
        #print 'rr', r

        
if __name__ == "__main__":
    unittest.main()

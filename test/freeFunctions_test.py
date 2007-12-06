#!/usr/bin/env python

__author__    = 'Koichi Takahashi <shafi@e-cell.org>'
__license__   = 'GPL'
__copyright__ = 'Copyright The Molecular Sciences Institute 2006-2007'


import unittest

import _gfrd as mod

import numpy


class FreeFunctionsTestCase( unittest.TestCase ):

    def setUp( self ):
        pass

    def tearDown( self ):
        pass


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
        result = scipy.integrate.quad( mod.g_bd, sigma, sigma + 
                                       6 * math.sqrt( 6 * D * t ),
                                       args=( sigma, t, D ) )
        igbd = result[0]
        print igbd
        self.failIf( ibd == 0 )
        self.assertAlmostEqual( 0.0, (ibd-igbd)/ibd )


    def test_int_g_bd_is_I_bd_smallt( self ):

        import scipy.integrate
        import math

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

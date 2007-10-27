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


#!/usr/bin/env python

__author__    = 'Koichi Takahashi <shafi@e-cell.org>'
__license__   = 'GPL'
__copyright__ = 'Copyright The Molecular Sciences Institute 2006-2007'


import unittest

import _gfrd as mod

import numpy


class FreeGreensFunctionTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_instantiation(self):
        D = 1e-12

        gf = mod.FreeGreensFunction(D)
        self.failIf(gf == None)

    def test_drawR(self):
        D = 1e-12
        
        gf = mod.FreeGreensFunction(D)

        t = 1e-3

        r = gf.drawR(0.5, t)
        self.failIf(r < 0.0)

        r = gf.drawR(0.0, t)
        self.failIf(r < 0.0)

        r = gf.drawR(1.0, t)
        self.failIf(r < 0.0)

    def test_drawR_zerot_is_zero(self):

        D = 1e-12
        
        gf = mod.FreeGreensFunction(D)

        t = 0.0

        r = gf.drawR(0.5, t)
        self.assertEqual(0.0, r)

        r = gf.drawR(0.0, t)
        self.assertEqual(0.0, r)

        r = gf.drawR(1.0, t)
        self.assertEqual(0.0, r)


    def no_test_ip_r_infinity_is_one(self):

        D = 1e-12
        
        t = 1e-5
        r = 2.5e-8
        
        gf = mod.FreeGreensFunction(D)
        ip = gf.ip_r(numpy.inf, t)
        self.assertEqual(1.0, ip)

    def test_int_p_r_is_ip_r(self):

        import scipy.integrate

        D = 1e-12
        t = 1e-5
        
        gf = mod.FreeGreensFunction(D)

        ip = gf.ip_r(0.0, t)
        self.assertEqual(0.0, ip)

        maxr = 5e-8

        resolution = 20
        for i in range(1, resolution):
            r = i * maxr / resolution 
            ip = gf.ip_r(r, t)
            result = scipy.integrate.quad(gf.p_r, 0.0, r,
                                          args=(t, ))
            np = result[0]
            self.assertAlmostEqual(0.0, (np-ip)/ip)


if __name__ == "__main__":
    unittest.main()

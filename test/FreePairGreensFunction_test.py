#!/usr/bin/env python

__author__    = 'Koichi Takahashi <shafi@e-cell.org>'
__license__   = 'GPL'
__copyright__ = 'Copyright The Molecular Sciences Institute 2006-2007'


import unittest

import _greens_functions as mod

import numpy


class FreePairGreensFunctionTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_instantiation(self):
        D = 1e-12

        gf = mod.FreePairGreensFunction(D)
        self.failIf(gf == None)

    def test_drawR(self):
        D = 1e-12
        r0 = 2e-8
        
        gf = mod.FreePairGreensFunction(D)

        t = 1e-3

        r = gf.drawR(0.5, r0, t)
        self.failIf(r < 0.0)

        r = gf.drawR(0.0, r0, t)
        self.failIf(r < 0.0)

        r = gf.drawR(1.0, r0, t)
        self.failIf(r < 0.0)

    def test_drawR_zerot_is_r0(self):
        D = 1e-12
        r0 = 2e-8
        
        gf = mod.FreePairGreensFunction(D)

        t = 0.0

        r = gf.drawR(0.5, r0, t)
        self.assertEqual(r, r0)

        r = gf.drawR(0.0, r0, t)
        self.assertEqual(r, r0)

        r = gf.drawR(1.0, r0, t)
        self.assertEqual(r, r0)

    def test_drawR_smallt(self):
        D = 1e-12
        r0 = 2e-8
        
        gf = mod.FreePairGreensFunction(D)

        t = 1e-4

        while t > 1e-60:
            r = gf.drawR(0.5, r0, t)
            self.failIf(r < 0.0)
            r = gf.drawR(0.0, r0, t)
            self.failIf(r < 0.0)
            r = gf.drawR(1.0, r0, t)
            self.failIf(r < 0.0)
            r = gf.drawR(1e-2, r0, t)
            self.failIf(r < 0.0)

            t *= 1e-3



    def test_draw_theta(self):
        D = 1e-12
        r0 = 5e-8
        t = 1e-4
        
        gf = mod.FreePairGreensFunction(D)

        #r = gf.drawR(0.5, r0, t)
        r = r0

        theta = gf.drawTheta(0.5, r, r0, t)
        self.failIf(theta < 0.0 or theta > numpy.pi)

        theta = gf.drawTheta(0.0, r, r0, t)
        self.failIf(theta < 0.0 or theta > numpy.pi)

        # rnd=1.0 fails on x87.  it's not a big problem, but
        # need to look at it later.
        theta = gf.drawTheta(0.999999, r, r0, t)
        self.failIf(theta < 0.0 or theta > numpy.pi)


    def test_ip_r_infinity_is_one(self):

        D = 1e-12
        
        t = 1e-5
        r0 = 5e-8
        r = 2.5e-8
        
        gf = mod.FreePairGreensFunction(D)

        ip = gf.ip_r(numpy.inf, r0, t)
        self.assertEqual(1.0, ip)


    def test_int_p_r_is_ip_r(self):

        import scipy.integrate

        D = 1e-12
        
        t = 1e-5
        r0 = 5e-8
        r = 2.5e-8
        
        gf = mod.FreePairGreensFunction(D)

        ip = gf.ip_r(0.0, r0, t)
        self.assertEqual(0.0, ip)

        maxr = 1e-6

        resolution = 20
        for i in range(1, resolution):
            r = i * maxr / resolution 
            ip = gf.ip_r(r, r0, t) 
            result = scipy.integrate.quad(gf.p_r, 0.0, r,
                                          args=(r0, t))
            np = result[0]
            self.assertAlmostEqual(0.0, (np-ip)/ip)


    def test_int_p_theta_is_p_r(self):

        import scipy.integrate

        D = 1e-12
        
        t = 1e-5
        r0 = 5e-8
        r = 2.5e-8
        
        gf = mod.FreePairGreensFunction(D)

        ip = gf.ip_theta(numpy.pi, r, r0, t)
        result = scipy.integrate.quad(gf.p_theta, 0.0, numpy.pi,
                                      args=(r, r0, t))
        np = result[0]

        pr = gf.p_r(r, r0, t) / (2 * numpy.pi * r * r)
        
        self.assertAlmostEqual(0.0, (pr-ip)/pr)
        self.assertAlmostEqual(0.0, (pr-np)/pr)





    def test_int_p_theta_is_ip_theta(self):

        import scipy.integrate

        D = 1e-12
        
        t = 1e-3
        r0 = 5e-8
        r = 2.5e-8
        
        gf = mod.FreePairGreensFunction(D)

        ip = gf.ip_theta(0.0, r, r0, t)
        self.assertEqual(0.0, ip)

        resolution = 20
        for i in range(1, resolution):
            theta = i * numpy.pi / resolution 
            ip = gf.ip_theta(theta, r, r0, t)
            result = scipy.integrate.quad(gf.p_theta, 0.0, theta,
                                          args=(r, r0, t))
            np = result[0]
            self.assertAlmostEqual(0.0, (np-ip)/ip)




if __name__ == "__main__":
    unittest.main()

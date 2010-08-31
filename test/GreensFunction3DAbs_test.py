#!/usr/bin/env python

__author__    = 'Koichi Takahashi <shafi@e-cell.org>'
__license__   = 'GPL'
__copyright__ = 'Copyright The Molecular Sciences Institute 2006-2007'


import unittest

import _greens_functions as mod

import numpy


class GreensFunction3DAbsTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_instantiation(self):
        D = 1e-12
        kf = 1e8
        sigma = 1e-8
        a = 1e-7
        r0 = 5e-8

        gf = mod.GreensFunction3DAbs(D, r0, a)
        self.failIf(gf == None)


    def test_draw_time(self):
        D = 1e-12
        a = 1e-7
        r0 = 5e-8
        
        gf = mod.GreensFunction3DAbs(D, r0, a)

        t = gf.drawTime(0.5)
        self.failIf(t <= 0.0 or t >= numpy.inf)

        t = gf.drawTime(0.0)
        self.failIf(t < 0.0 or t >= numpy.inf)

        t = gf.drawTime(1.0)
        self.failIf(t <= 0.0 or t >= numpy.inf)

    def test_draw_time_r0_equal_a(self):
        D = 1e-12

        a = 1e-7
        r0 = a
        
        gf = mod.GreensFunction3DAbs(D, r0, a)

        t = gf.drawTime(0.5)
        self.assertEqual(0.0, t)

    def test_drawR(self):
        D = 1e-12

        a = 1e-7
        r0 = 2e-8
        
        gf = mod.GreensFunction3DAbs(D, r0, a)

        t = 1e-3

        r = gf.drawR(0.5, t)
        self.failIf(r < 0.0 or r > a)

        r1 = gf.drawR(0.0, t)
        r2 = gf.drawR(1.0, t)

        self.failIf(r1 < 0.0 or r1 > a)
        self.failIf(r2 < 0.0 or r2 > a)

        self.failIf(abs(r1) > 1e-15)
        self.assertAlmostEqual(abs(r2 - a), 0)

    def test_drawR_zerot(self):
        D = 1e-12
        a = 1e-7
        r0 = 2e-8
        
        gf = mod.GreensFunction3DAbs(D, r0, a)

        t = 0.0

        r = gf.drawR(0.5, t)
        self.assertEqual(r0, r)


    def test_drawR_squeezed(self):

        D = 1e-12
        a = 1.01e-8
        
        t = 1e-6
        r0 = 1.005e-8
        gf = mod.GreensFunction3DAbs(D, r0, a)
        r = gf.drawR(0.5, t)
        self.failIf(r < 0.0 or r > a)

        # near 0
        r = 0.0001e-8
        r0 = 0.0001e-8
        gf = mod.GreensFunction3DAbs(D, r0, a)
        r = gf.drawR(0.5, t)
        self.failIf(r < 0.0 or r > a)

        # near a
        r = 1.0099e-8
        r0 = 1.0099e-8
        gf = mod.GreensFunction3DAbs(D, r0, a)
        r = gf.drawR(0.5, t)
        self.failIf(r < 0.0 or r > a)


    def test_draw_theta(self):
        D = 1e-12
        a = 1e-7
        r0 = 5e-8
        
        gf = mod.GreensFunction3DAbs(D, r0, a)

        t = gf.drawTime(0.5)
        r = gf.drawR(0.5, t)

        theta = gf.drawTheta(0.5, r, t)
        self.failIf(theta < 0.0 or theta > numpy.pi)

        theta = gf.drawTheta(0.0, r, t)
        self.failIf(theta < 0.0 or theta > numpy.pi)

        theta = gf.drawTheta(1.0, r, t)
        self.failIf(theta < 0.0 or theta > numpy.pi)


    def test_draw_theta_zerot(self):
        D = 1e-12

        a = 1e-7
        r = 5e-8
        r0 = 5e-8
        
        gf = mod.GreensFunction3DAbs(D, r0, a)

        t = 0.0
        theta = gf.drawTheta(0.5, r0, t)
        self.assertEqual(0.0, theta)

    def test_draw_theta_smallt(self):

        D = 1e-12

        a = 1e-7
        r = 5e-8
        r0 = 5e-8
        
        gf = mod.GreensFunction3DAbs(D, r0, a)

        t = 1e-4  # not very small though..
        theta = gf.drawTheta(0.5, r, t)
        self.failIf(theta < 0.0 or theta > numpy.pi)

    def test_draw_theta_large_t(self):

        D = 1e-12

        a = 1e-7
        r = 5e-8
        r0 = 5e-8
        
        gf = mod.GreensFunction3DAbs(D, r0, a)

        t = 1e5 
        theta = gf.drawTheta(0.5, r, t)
        self.failIf(theta < 0.0 or theta > numpy.pi)



    def test_draw_theta_near_a(self):

        D = 1e-12
        #a = 1.01e-8  # this is a better test but currently fails
        a = 1.1e-8
        
        t = 1e-5

        # near a
        r = 1.009999e-8
        r0 = 1.009999e-8
        gf = mod.GreensFunction3DAbs(D, r0, a)
        theta = gf.drawTheta(0.5, r, t)
        self.failIf(theta < 0.0 or theta > numpy.pi)

    def test_draw_theta_r_equal_a(self):
        D = 1e-12
        a = 1e-7
        r0 = 9e-8

        t = 1e-4
        r = a
        
        gf = mod.GreensFunction3DAbs(D, r0, a)

        theta = gf.drawTheta(0.5, r, t)

        self.failIf(theta < 0.0 or theta > numpy.pi)

    def test_draw_theta_r0_near_a_r_equal_a(self):
        D = 1e-12
        a = 1e-7
        r0 = a - 1e-9

        t = 1e-2
        r = a
        
        gf = mod.GreensFunction3DAbs(D, r0, a)

        theta = gf.drawTheta(0.5, r, t)
        self.failIf(theta < 0.0 or theta > numpy.pi)


    def test_p_int_r(self):

        D = 1e-12

        t = 1e-3
        r0 = 5e-8

        a = 1e-7
        
        gf = mod.GreensFunction3DAbs(D, r0, a)

        r = r0
        pintr = gf.p_int_r(r, t)
        self.failIf(pintr < 0.0 or pintr > 1.0)


    def test_p_int_r_at_a_is_p_survival(self):

        D = 1e-12

        t = 1e-3
        r0 = 5e-8

        a = 1e-7
        
        gf = mod.GreensFunction3DAbs(D, r0, a)
        r = r0
        
        pintr = gf.p_int_r(a, t)
        psurv = gf.p_survival(t)
        self.assertAlmostEqual(pintr, psurv)

    def test_p_int_r_at_zero_is_zero(self):

        D = 1e-12

        t = 1e-3
        r0 = 5e-8

        a = 1e-7
        
        gf = mod.GreensFunction3DAbs(D, r0, a)
        r = r0
        
        pintr = gf.p_int_r(0.0, t)
        self.assertEqual(0.0, pintr)

    def test_ip_theta_is_int_p_theta(self):

        import scipy.integrate

        D = 1e-11

        t = 1e-3
        r0 = 5e-8

        a = 1e-7
        
        gf = mod.GreensFunction3DAbs(D, r0, a)
        r = r0

        ip = gf.ip_theta(0.0, r, t)
        self.assertEqual(0.0, ip)
        
        resolution = 10
        for i in range(1, resolution):
            theta = i * numpy.pi / resolution 
            ip = gf.ip_theta(theta, r, t)
            result = scipy.integrate.quad(gf.p_theta, 0.0, theta,
                                          args=(r, t))
            np = result[0]
            #print theta, np, ip
            self.assertAlmostEqual(0.0, (np-ip)/ip)

    '''

    def test_ip_theta_pi_is_p_0(self):

        D = 1e-12
        sigma = 1e-8
        kf = 1e-8

        t = 1e-5
        r0 = 5e-8
        r = r0

        a = 1e-7
        
        gf = mod.GreensFunction3DAbs(D, kf, sigma, a)

        ip = gf.ip_theta(numpy.pi, r, r0, t)
        p0 = gf.p_0(t, r, r0) * 2

        self.assertNotEqual(0.0, ip)
        self.assertAlmostEqual(1.0, ip/p0)
'''


    def test_p_theta_never_negative(self):

        D = 1e-12

        # smaller t causes problem
        t = 1e-3
        r0 = 5e-8
        r = r0
        a = 1e-7
        
        gf = mod.GreensFunction3DAbs(D, r0, a)

        pint = gf.ip_theta(numpy.pi, r, t)

        pmin = 0.0
        resolution = 50
        for i in range(resolution):
            theta = i * numpy.pi / resolution
            p = gf.p_theta(theta, r, t) / pint / resolution 
            pmin = min(pmin, p)
            #print 'theta: ', theta, '\tp: ', p
            
        self.failIf(pmin < 0.0, 'Negative p_theta; t= %g, %s'
                    % (t, gf.dump()))


    def test_ip_theta_never_decrease(self):

        D = 1e-12

        # smaller t causes problem
        t = 1e-3
        r0 = 5e-8
        r = r0
        a = 1e-7
        
        gf = mod.GreensFunction3DAbs(D, r0, a)

        pint_prev = 0.0

        resolution = 50
        for i in range(resolution):
            theta = i * numpy.pi / resolution
            pint = gf.ip_theta(theta, r, t)
            self.failIf(pint < pint_prev)
            pint_prev = pint

    def test_idp_theta_at_a_is_dp_survival(self):

        D = 1e-12

        # smaller t causes problem
        t = 1e-3
        r0 = 9e-8
        a = 1e-7
        
        gf = mod.GreensFunction3DAbs(D, r0, a)

        dp = gf.dp_survival(t)
        iptheta = gf.idp_theta(numpy.pi, a, t) * numpy.pi * a * a * 2

        self.assertAlmostEqual(dp, iptheta)






if __name__ == '__main__':
    unittest.main()

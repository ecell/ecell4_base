#!/usr/bin/env python

import unittest

import numpy

import _greens_functions as mod

class GreensFunction3DAbsSymTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_instantiation(self):
        D = 1e-12
        a = 1.0
        gf = mod.GreensFunction3DAbsSym(D, a)
        self.failIf(gf == None)

    def test_no_shell(self):
        D = 1e-12
        a = numpy.inf
        gf = mod.GreensFunction3DAbsSym(D, a)
        t = gf.drawTime(0.5)
        self.assertEqual(numpy.inf, t)

        # not supported now
        #r = gf.drawR(0.5, 1e-6)
        #self.assertAlmostEqual(p_free, r)


    def test_zero_shell(self):
        D = 1e-12
        a = 0.0
        gf = mod.GreensFunction3DAbsSym(D, a)
        
        t = gf.drawTime(0.5)
        self.assertEqual(0.0, t)
        r = gf.drawR(0.5, t)
        self.assertEqual(0.0, r)


    def test_draw_time(self):
        D = 1e-12
        a = 1e-7
        gf = mod.GreensFunction3DAbsSym(D, a)

        t = gf.drawTime(0.0)
        t = gf.drawTime(0.5)
        t = gf.drawTime(1.0 - 1e-16)

    def test_drawR(self):
        D = 1e-12
        a = 1e-7
        gf = mod.GreensFunction3DAbsSym(D, a)

        t = gf.drawTime(0.5)

        r = gf.drawR(0.0, t)
        r = gf.drawR(0.5, t)
        r = gf.drawR(1.0 - 1e-16, t)

    def test_drawR_zerot(self):
        D = 1e-12
        a = 1e-8
        t = 0.0

        gf = mod.GreensFunction3DAbsSym(D, a)

        r = gf.drawR(0.5, t)
        self.assertEqual(0.0, r)

    def test_drawR_smallt(self):
        D = 1e-12
        a = 1e-8
        gf = mod.GreensFunction3DAbsSym(D, a)
        t = gf.drawTime(0.5)

        while t > 1e-30:
            t *= 1e-4
            r = gf.drawR(0.0, t)
            self.failIf(r < 0.0)
            self.failIf(r > a)


    def test_drawR_large_shell(self):
        D = 1e-12
        a = 1e-3
        gf = mod.GreensFunction3DAbsSym(D, a)
        t = 1e-10
        r = gf.drawR(0.5, t)
        self.failIf(r <= 0.0)
        self.failIf(r > a)

    def test_drawR_large_t(self):
        D = 1e-12
        a = 1e-6
        gf = mod.GreensFunction3DAbsSym(D, a)
        t = gf.drawTime(0.0)
        r = gf.drawR(0.5, t)

        self.failIf(r <= 0.0)
        self.failIf(r > a)


    def test_drawR_regression1(self):

        # this parameter set causes very small p_survival.

        D = 1e-12
        a=3.30588e-08
        rnd=0.944502 
        dt=0.00411832

        gf = mod.GreensFunction3DAbsSym(D, a)
        r = gf.drawR(rnd, dt)

        self.failIf(r <= 0.0)
        self.failIf(r > a)


    def test_p_int_r_is_p_int_r_free_with_large_shell(self):
        D = 1e-12
        a = 1e-6
        gf = mod.GreensFunction3DAbsSym(D, a)

        r = 1e-9
        t = 1e-6

        p = gf.p_int_r(r, t)
        p_free = gf.p_int_r_free(r, t)

        self.assertAlmostEqual(p, p_free)


    def test_p_int_r_at_a_is_p_survival(self):
        D = 1e-12
        a = 1e-8
        gf = mod.GreensFunction3DAbsSym(D, a)

        t = 1e-5

        p = gf.p_int_r(a, t)
        psurv = gf.p_survival(t)
        print 'ps', p, psurv
        self.assertAlmostEqual(p, psurv)

    def test_p_int_r_at_a_is_p_survival_2(self):
        D = 1e-12
        a = 1e-9
        gf = mod.GreensFunction3DAbsSym(D, a)

        t = 1e-3

        p = gf.p_int_r(a, t)
        psurv = gf.p_survival(t)
        print 'ps', p, psurv
        self.assertAlmostEqual(p, psurv)




if __name__ == "__main__":
    unittest.main()

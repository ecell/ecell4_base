#!/usr/bin/env python

import math
import numpy

import _gfrd



def p_rev(r, t, r0, kf, D, sigma):

    jacobian = 4.0 * numpy.pi * r * r

    








'''
x + y + z - h == 0
x * y + y * z + x * z - kd == 0
x * y * z - a == 0
'''

class Func:

    def __init__(self, h, kd, a):
        self.h = h
        self.kd = kd
        self.a = a

    def __call__(self, x):
        h = self.h
        kd = self.kd
        a = self.a

        result = numpy.array([x[0] + x[1] + x[2] - h,
                              x[0] * x[1] + x[1] * x[2] + x[0] * x[2] - kd,
                              x[0] * x[1] * x[2] - a])
        return result

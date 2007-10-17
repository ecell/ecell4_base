#!/usr/bin/env python

import math
import numpy

import _gfrd



def p_irr( r, t, r0, kf, D, sigma ):

    return _gfrd.p_irr_radial( r, t, r0, kf, D, sigma )




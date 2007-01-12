

import gfrd
import math
import random

import numpy


t = 1e-4
D = 1e-12
a = 2e-8


gf = gfrd.FirstPassageGreensFunction( D )

#print gf.p_r_fourier( 1e-10, 1e-6 )

print gf.p_r_int( 2e-9, t, a )
print gf.p_survival( t, a )

#print gf.drawTime( .5, a )
#for i in range(100000):
#    gf.drawTime( random.random(), a )

print gf.drawR( .5, t, a )
for i in range(100000):
    gf.drawR( random.random(), t, a )

#print gf.p_r_fourier( 1e-9, 1e-6 )
#for i in range(100000):
#    gf.p_r_fourier( 1e-10, 1e-6 )

#for i in range(1000):
#    numpy.random.normal( 0.0, 1, 3 )

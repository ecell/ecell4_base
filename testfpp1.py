

import gfrd
import math

import numpy


t = 1e-4
D = 1e-12
a = 2e-9


gf = gfrd.FirstPassageGreensFunction( D, a )

#print gf.p_r_fourier( 1e-10, 1e-6 )

print gf.p_r_int( 2e-9, 5e3 )
print gf.p_survival( 5e3 )

#print gf.p_r_fourier( 1e-9, 1e-6 )
#for i in range(100000):
#    gf.p_r_fourier( 1e-10, 1e-6 )

#for i in range(1000):
#    numpy.random.normal( 0.0, 1, 3 )

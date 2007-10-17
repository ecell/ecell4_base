#!/usr/bin/env/python


import numpy
from matplotlib.pylab import *

from p_irr import p_irr

N = 100

sigma = 1e-9
r0 = sigma
t = 1e-5
kf = 1e-10
D = 1e-11

rmin = sigma
rmax = sigma * 1e2

rtick = ( rmax - rmin ) / N
rarray = numpy.array( range( N ) ) * rtick + rmin
parray = numpy.zeros( N )

for i in range( N ):
    #print p_irr( rarray[i], t, r0, kf, D, sigma )
    parray[i] = p_irr( rarray[i], t, r0, kf, D, sigma )
    
#print rarray, parray

plot( parray )
show()


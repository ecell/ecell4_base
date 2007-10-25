#!/usr/bin/env/python

import sys

import numpy
import scipy.io
from matplotlib.pylab import *


import _gfrd
from p_irr import p_irr

infilename = sys.argv[1]


N_A = 6.0221367e23

N = 100

sigma = 1e-7
r0 = sigma + 1e-18
t = .1
kf = 1e6 / N_A
D = 1e-11

rmin = sigma
rmax = sigma * 1e2

rtick = ( rmax - rmin ) / N
rarray = numpy.array( range( N ) ) * rtick + rmin
parray = numpy.zeros( N )


for i in range( N ):
    parray[i] = p_irr( rarray[i], t, r0, kf, D, sigma )

S = _gfrd.S_irr( t, r0, kf, D, sigma )

parray /= parray.sum()
parray *= S 

print rarray, parray

plot( rarray, parray  )



bins = 100

data = scipy.io.read_array( infilename )

#grid = numpy.mgrid[rmin:rmax:rmax/N2]
nonreactions = numpy.compress( data >= sigma, data )
hist, lower_edges = numpy.histogram( nonreactions, bins=bins )

histsum = hist.sum()
S_sim = float( len( nonreactions ) ) / len( data )
hist = hist.astype( numpy.float )
hist /= len( data )


x = lower_edges + (lower_edges[1]-lower_edges[0]) * .5

plot( lower_edges, hist )

print S, S_sim


show()


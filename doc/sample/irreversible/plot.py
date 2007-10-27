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
rarray = numpy.mgrid[rmin:rmax:rtick]
parray = numpy.zeros( N )


for i in range( N ):
    parray[i] = p_irr( rarray[i], t, r0, kf, D, sigma )

S = _gfrd.S_irr( t, r0, kf, D, sigma )

#parray /= parray.sum() * rtick
parray *= S 

print rarray, parray

bins = 100
print 'load'
#data = scipy.io.read_array( infilename )  # <-- slow
infile = open( infilename )
data = array([float(x) for x in infile.read().split()], numpy.float)
infile.close()

print 'hist'
#grid = numpy.mgrid[rmin:rmax:rmax/N2]
nonreactions = numpy.compress( data >= sigma, data )
hist, lower_edges = numpy.histogram( nonreactions, bins=bins )

histsum = hist.sum()
S_sim = float( len( nonreactions ) ) / len( data )
hist = hist.astype( numpy.float )

print len( lower_edges[1:] ), len( lower_edges[:-1])
print len( hist )
xtick = lower_edges[2]-lower_edges[1]

hist /= len( data )
hist /= xtick
#hist *= 3.14

x = lower_edges + ( xtick * .5 )

print hist,lower_edges, lower_edges[1:] - lower_edges[:-1]
plot( rarray, parray  )
plot( x, hist )

print S, S_sim


show()


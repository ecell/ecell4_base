#!/usr/bin/env/python

import sys

import numpy
import scipy.io
from matplotlib.pylab import *


import _gfrd
from p_rev import p_rev

infilename = sys.argv[1]


N_A = 6.0221367e23

N = 100

sigma = 1e-7
r0 = sigma
t = .01
kf = 1e6 / N_A
D = 1e-11

# rmin = sigma
# rmax = sigma * 1e2
# rtick = ( rmax - rmin ) / N
# rarray = numpy.mgrid[rmin:rmax:rtick]
# parray = array( [ p_rev( r, t, r0, kf, D, sigma ) for r in rarray ] )
# S = _gfrd.S_irr( t, r0, kf, D, sigma )

data = scipy.io.read_array( 'p_rev.tsv' )
print data
rarray, parray = numpy.transpose( data )

print rarray, parray

bins = 50
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

hist /= len( data ) * xtick
#hist /= xtick
#hist *= 3.14

x = lower_edges + ( xtick * .5 )

print hist,lower_edges, lower_edges[1:] - lower_edges[:-1]
plot( rarray, parray  )
plot( x, hist )

#print S, S_sim


show()


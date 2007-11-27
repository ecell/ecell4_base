#!/usr/bin/env/python

import sys

import numpy
import scipy.io
from matplotlib.pylab import *

#import _gfrd

infilename = sys.argv[1]


N_A = 6.0221367e23

N = 100

sigma = 1e-7

r0 = sigma
kf = 1e6 / N_A
D = 1e-11

#tau = sigma*sigma / D
#t = .01


data = scipy.io.read_array( 'p_rev.tsv' )
rarray, parray = numpy.transpose( data )

rarray /= sigma

print rarray, parray

bins = 25
print 'load'
#data = scipy.io.read_array( infilename )  # <-- slow
infile = open( infilename )
data = array([float(x) for x in infile.read().split()], numpy.float)
infile.close()

print 'hist'
nonreactions = numpy.compress( data >= sigma, data )
hist, lower_edges = numpy.histogram( nonreactions, bins=bins )

histsum = hist.sum()
S_sim = float( len( nonreactions ) ) / len( data )
hist = hist.astype( numpy.float )

xtick = lower_edges[2]-lower_edges[1]

hist /= len( data ) * xtick

x = lower_edges + ( xtick * .5 )
x /= sigma

print hist,lower_edges, lower_edges[1:] - lower_edges[:-1]
plot( rarray, parray, 'b-'  )
plot( x, hist, 'k.' )
loglog( rarray, parray, 'b-'  )
loglog( x, hist, 'k.' )

show()


#!/usr/bin/env/python

import sys

import numpy
import scipy.io
from matplotlib.pylab import *


import gfrdbase

N_A = gfrdbase.N_A


D = 1e-11

sigma = 1e-7
r0 = sigma
kf = 1e6 / N_A
t = 0.01




def plot_sol( rmax ):
    rmin = sigma
    
    N = 100
    rtick = ( rmax - rmin ) / N
    rarray = numpy.mgrid[rmin:rmax:rtick]
    
    parray = array( [ gfrdbase.p_free( r, t, D ) for r in rarray ] )
    
    print rarray, parray

    plot( rarray / sigma , parray, 'b-'  )



def plot_file( infilename, maxr ):
    
    bins = 25
    print 'load'
    infile = open( infilename )
    data = array([float(x) for x in infile.read().split()], numpy.float)
    infile.close()

    print 'hist'
    nonreactions = numpy.compress( data >= sigma, data )
    print 'max', max( nonreactions )
    hist, lower_edges = numpy.histogram( nonreactions, bins=bins,
                                         range=[ sigma, maxr ] )
    print 'hist', hist
    print 'le', lower_edges
    
    histsum = hist.sum()
    S_sim = float( len( nonreactions ) ) / len( data )
    hist = hist.astype( numpy.float )
    
    xtick = lower_edges[2]-lower_edges[1]
    
    hist /= len( data ) * xtick
    
    x = lower_edges + ( xtick * .5 )
    print 'x', x

    plot( x / sigma, hist, '.', label=infilename )


if __name__ == '__main__':

    plot_sol( sigma * 25 )

    for filename in sys.argv[1:]:
        plot_file( filename, sigma * 25 )
    
    legend()
    show()



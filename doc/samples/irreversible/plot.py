#!/usr/bin/env/python

import sys

import numpy
import scipy.io
from matplotlib.pylab import *


import _gfrd
from p_irr import p_irr

N_A = 6.0221367e23

N = 100

sigma = 1e-7
r0 = sigma
kf = 1e6 / N_A
D = 1e-11

tau = sigma*sigma / D

rmin = sigma

def load_data( filename ):
    infile = open( filename )
    data = array([float(x) for x in infile.read().split()], numpy.float)
    infile.close()

    return data
    
def plot_sol( t, rmax ):

    rtick = ( rmax - rmin ) / N
    rarray = numpy.mgrid[rmin:rmax:rtick]

    parray = array( [ p_irr( r, t, r0, kf, D, sigma ) for r in rarray ] )

    loglog( rarray / sigma , parray, '-', label='theory' )

def plot_hist( data, T ):

    bins = 25

    nonreactions = numpy.compress( data >= sigma, data )
    print 'max', max( nonreactions )
    hist, lower_edges = numpy.histogram( nonreactions, bins=bins )

    histsum = hist.sum()
    S_sim = float( len( nonreactions ) ) / len( data )
    hist = hist.astype( numpy.float )

    xtick = lower_edges[2]-lower_edges[1]

    hist /= len( data ) * xtick

    x = lower_edges + ( xtick * .5 )
    #print 'x', x

    loglog( x / sigma, hist, '.', label='sim (T = %g tau)' % (T * 100) )
    


if __name__ == '__main__':

    for i in range( len(sys.argv[1:])/2 ):
        filename = sys.argv[i*2+1]
        T = float( sys.argv[i*2+2] )
        print filename,T
        data = load_data( filename )
        plot_sol( T, max( data ) )
        plot_hist( data, T )


    xlabel( 'r / sigma' )
    ylabel( 'p_irr' )
    legend()
    show()




#S = _gfrd.S_irr( t, r0, kf, D, sigma )
#print S, S_sim



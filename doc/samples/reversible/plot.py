#!/usr/bin/env/python

# PYTHONPATH=../../.. python plot.py rev.-1.out p_rev.-1.tsv 0.00001 rev.0.out p_rev.0.tsv 0.0001 rev.1.out p_rev.1.tsv 0.001 rev.2.out p_rev.2.tsv 0.01 
# rev.3.out p_rev.3.tsv 0.1

import sys

import numpy
import scipy.io
from matplotlib.pylab import *

#import _gfrd

infilename = sys.argv[1]


N_A = 6.0221367e23

N = 100

sigma = 1e-8

#r0 = sigma
#D = 1e-12
#kf = 10 * sigma * D

#tau = sigma*sigma / D
#t = .01


def load_data( filename ):
    infile = open( filename )
    data = array([float(x) for x in infile.read().split()], numpy.float)
    infile.close()

    return data
    
def plot_sol( filename, maxr ):

    data = scipy.io.read_array( filename )
    rarray, parray = numpy.transpose( data )
    mask = numpy.less_equal( rarray, maxr )
    rarray = numpy.compress( mask, rarray )
    parray = numpy.compress( mask, parray )
    print rarray, parray

    loglog( rarray / sigma, parray, '-'  )



def plot_hist( data, T ):

    bins = 20

    nonreactions = numpy.compress( data > sigma, data )
    print 'max', max( nonreactions )
    hist, lower_edges = numpy.histogram( nonreactions, bins=bins )

    histsum = hist.sum()
    S_sim = float( len( nonreactions ) ) / len( data )
    hist = hist.astype( numpy.float )
    print hist
    xtick = lower_edges[2]-lower_edges[1]

    hist /= len( data ) * xtick

    x = lower_edges + ( xtick * .5 )
    print 'x', x, hist

    loglog( x / sigma, hist, '.', label='sim (T = %g tau)' % (T * 100) )
    


if __name__ == '__main__':

    for i in range( len(sys.argv[1:])/3 ):
        simfilename = sys.argv[i*3+1]
        solfilename = sys.argv[i*3+2]
        T = float( sys.argv[i*3+3] )
        print simfilename,solfilename,T
        data = load_data( simfilename )
        plot_sol( solfilename, max( data ) )
        plot_hist( data, T )


    xlabel( 'r / sigma' )
    ylabel( 'p_rev' )
    #legend()
    show()



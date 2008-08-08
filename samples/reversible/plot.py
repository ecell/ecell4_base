#!/usr/bin/env/python

# PYTHONPATH=../.. python plot.py rev.-3.out p_rev.-3.tsv 0.0000000125 rev.-2.out p_rev.-2.tsv 0.000000125 rev.-1.out p_rev.-1.tsv 0.00000125 rev.0.out p_rev.0.tsv 0.0000125 rev.1.out p_rev.1.tsv 0.000125 rev.2.out p_rev.2.tsv 0.00125 rev.3.out p_rev.3.tsv 0.0125

import sys

import numpy
import scipy.io
from matplotlib.pylab import *

#import _gfrd

infilename = sys.argv[1]


N_A = 6.0221367e23

sigma = 5e-9

#r0 = sigma
D_tot = 2e-12
#kf = 10 * sigma * D

tau = sigma*sigma / D_tot
rmin = sigma


def load_data( filename ):
    infile = open( filename )
    data = array([float(x) for x in infile.read().split()], numpy.float)
    infile.close()

    return data
    
def plot_sol( filename, t ):

    rmax = 3.1 * math.sqrt( 6 * D_tot * t ) + rmin

    data = scipy.io.read_array( filename )
    rarray, parray = numpy.transpose( data )
    mask = numpy.less_equal( rarray, rmax )
    rarray = numpy.compress( mask, rarray )
    parray = numpy.compress( mask, parray )

    return loglog( rarray / sigma, parray, 'k-'  )[0]



def plot_hist( data, T, i ):

    bins = 20

    nonreactions = numpy.compress( data > sigma, data )
    print 'max', max( nonreactions )
    hist, lower_edges = numpy.histogram( nonreactions, bins=bins )

    histsum = hist.sum()
    S_sim = float( len( nonreactions ) ) / len( data )
    hist = hist.astype( numpy.float )
    #print hist
    xtick = lower_edges[2]-lower_edges[1]

    hist /= len( data ) * xtick

    x = lower_edges + ( xtick * .5 )
    #print 'x', x, hist
    colors = [ 'b', 'g', 'r', 'c', 'm', 'y', 'k' ]

    loglog( x / sigma, hist, colors[i] + '.',             
            label=r'$T = \tau^{%d}$' % round(math.log10(T/tau)) )

    return lower_edges[-1] + xtick
    


if __name__ == '__main__':

    for i in range( len(sys.argv[1:])/3 ):
        simfilename = sys.argv[i*3+1]
        solfilename = sys.argv[i*3+2]
        T = float( sys.argv[i*3+3] )
        print simfilename,solfilename,T
        data = load_data( simfilename )
        plot_hist( data, T, i )
        solline = plot_sol( solfilename, T )



    xlabel( r'$r / \sigma$', fontsize='large' )
    ylabel( r'$p_{rev}$', fontsize='large' )
    xlim( 0.9, 5e2 )
    ylim( 1.5e1, 7e9 )
    solline.set_label( r'theory' )
    legend( handlelen=0.02, pad=0.02,handletextsep=0.01, labelsep=0.001 )
    grid()
    show()


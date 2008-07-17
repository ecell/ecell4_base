#!/usr/bin/env python

#  PYTHONPATH=../.. python plot.py irr.-3.out 0.0000000125 irr.-2.out 0.000000125  irr.-1.out 0.00000125  irr.0.out 0.0000125 irr.1.out 0.000125 irr.2.out 0.00125 irr.3.out 0.0125


import sys

import numpy
import scipy.io
from matplotlib.pylab import *


import _gfrd
from p_irr import p_irr

N_A = 6.0221367e23

N = 10000

sigma = 5e-9
r0 = sigma
D_tot = 2e-12
kf = 10 * sigma *  D_tot
#kf = 0

tau = sigma*sigma / D_tot

rmin = sigma

def load_data( filename ):
    infile = open( filename )
    data = array([float(x) for x in infile.read().split()], numpy.float)
    infile.close()

    return data
    
def plot_sol( t ):

    rmax = 3.1 * math.sqrt( 6 * D_tot * t ) + rmin

    rtick = ( rmax - rmin ) / N
    rarray = numpy.mgrid[rmin:rmax:rtick]

    parray = array( [ p_irr( r, t, r0, kf, D_tot, sigma ) for r in rarray ] )

    return loglog( rarray / sigma , parray, 'k-' )[0]
    #plot( rarray / sigma , parray, 'k-', label='theory' )

def plot_hist( data, T, i ):

    bins = 20

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
    #pStyles = [ 'o', '^', 'v', '<', '>', 's', '+' ]
    colors = [ 'b', 'g', 'r', 'c', 'm', 'y', 'k' ]

    loglog( x / sigma, hist, colors[i] + '.', 
            label=r'$T = \tau^{%d}$' % round(math.log10(T/tau)) )

    
    return lower_edges[-1] + xtick


if __name__ == '__main__':

    for i in range( len(sys.argv[1:])/2 ):
        filename = sys.argv[i*2+1]
        T = float( sys.argv[i*2+2] )
        print filename,T
        data = load_data( filename )
        maxr = plot_hist( data, T, i )
        solline = plot_sol( T )



    xlabel( r'$r / \sigma$', fontsize='large' )
    ylabel( r'$p_{irr}$', fontsize='large' )
    xlim( 0.9, 5e2 )
    ylim( 1.5e1, 7e9 )
    solline.set_label( r'theory' )
    legend( handlelen=0.02, pad=0.02,handletextsep=0.01, labelsep=0.001 )
    grid()
    show()




#>>> _gfrd.S_irr( .0001 * 1e-8**2/1e-12, 1e-8, 10 * 1e-8 * 1e-12, 1e-12, 1e-8 )
#0.99116163945434221



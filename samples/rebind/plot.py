#!/usr/bin/env/python

# PYTHONPATH=../.. python plot.py rev.-3.out p_rev.-3.tsv 0.0000000125 rev.-2.out p_rev.-2.tsv 0.000000125 rev.-1.out p_rev.-1.tsv 0.00000125 rev.0.out p_rev.0.tsv 0.0000125 rev.1.out p_rev.1.tsv 0.000125 rev.2.out p_rev.2.tsv 0.00125 rev.3.out p_rev.3.tsv 0.0125

import sys

import numpy
import scipy.io
from matplotlib.pylab import *

#import _gfrd

infilename = sys.argv[1]


N_A = 6.0221367e23

N = 100

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
    

def plot_hist( data, xmin, xmax, N ):

    #    xmin = data.min()
    #xmax = data.max()

    logxmin = math.log(xmin)
    logxmax = math.log(xmax)
    
    tick=(logxmax-logxmin)/N
    loggrid = numpy.mgrid[logxmin:logxmax:tick]
    grid = numpy.exp(loggrid)
    print len(data)
    print grid, xmin,xmax
    
    n, bins = numpy.histogram(numpy.log10(data), bins=N)
    n = n.astype(numpy.floating)
    n /= float(len(data))
    x, y = 10**bins[:-1], n+1e-10
    print x.shape, y.shape
    loglog( x, y )#, label=filename )


if __name__ == '__main__':

    import os
    import glob

    xmin = 1e-9
    xmax = 1e3

    for i in range( len(sys.argv[1:])/1 ):
        simpattern = sys.argv[i+1]

        globpattern = simpattern.replace('ALL','*')
        l = os.path.basename( os.path.splitext( simpattern )[0] )
        print 'pattern ', l
        filelist = glob.glob( globpattern )
        print filelist
        data = []
        for file in filelist:
            data.append( load_data( file ) )
        data = numpy.concatenate( data )
        print len(data)

        plot_hist( data, xmin, xmax, N )


    xticks( [1e-15,1e-12, 1e-9, 1e-6, 1e-3, 1, 1e3], 
            [r'${\rm 1 fs}$',
             r'${\rm 1 ps}$',
             r'${\rm 1 ns}$',
             r'${\rm 1 \mu s}$',
             r'${\rm 1 ms}$',
             r'${\rm 1 s}$',
             r'${\rm 1000 s}$'],
            size=24 )
            #yticks( [],[] )

    #xlabel( r'$r / \sigma$', fontsize='large' )
    #xlabel( r'time [s]', fontsize='large' )
    ylabel( r'relative frequency', size=24 )
    #xlim( 1e-15, 10 )
    ylim( 2e-6, 1 )
    #solline.set_label( r'theory' )
    legend( handlelen=0.02, pad=0.02,handletextsep=0.01, labelsep=0.001 )
    grid()
    show()


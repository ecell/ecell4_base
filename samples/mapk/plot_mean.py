#!/usr/bin/env python

import sys

import numpy
import scipy.io

from matplotlib.pylab import *

def load_header( filename ):
    file = open( filename )
    header = []
    for line in file.readlines():
        if line[0:2] == '#@':
            hline = line[2:].lstrip()
            header.append( hline )

    return header

def resample( x, y, newx ):

    indices = numpy.searchsorted( x, newx )

    return y.take( indices )
    

def add_columns( data, ycolumns ):

    y = numpy.array([ data[:,col] for col in ycolumns ]) 

    y = y.sum(0)

    return y

def load_data( filename ):
    ycolumns = [1,]
    #ycolumns = [2,6]
    #ycolumns = [3,5]
    #ycolumns = [2,6,3,5]

    header = load_header( filename )
    print header
    for l in header:
        exec( l )

    data = numpy.loadtxt( filename )
    x = data[:,0]
    y = add_columns( data, ycolumns )

    return x, y


def plot_file( filename ):

    x, y = load_data( filename )

    #plot_theory( N_K, N_P, Keq, x[-1] )
    plot( x, y, '-' )

    #psd( y )
    #ylim( 1, 5e4 )


import glob
import os

for pattern in sys.argv[1:]:

    globpattern = pattern.replace('ALL','*')

    l = os.path.basename( os.path.splitext( pattern )[0] )
    print 'pattern ', l

    filelist = glob.glob( globpattern )

    start = 0
    end = 80.
    interval = (end-start) / 1000
    rx = numpy.mgrid[start:end:interval]

    data = []

    for filename in filelist:
        print 'file ', filename
        x, y = load_data( filename )
        ry = resample( x, y, rx )
        print ry.shape
        data.append( ry )

        mry = numpy.array( data ).mean( 0 )


    plot( rx, mry, label=l )

from matplotlib.font_manager import FontProperties
legend( loc='lower right', prop=FontProperties( size='tiny' ),pad=0.01 )

plot_file('/home/shafi/wrk/brown/samples/mapk/Kpp.ecd' )
plot_file('/home/shafi/wrk/brown/samples/mapk/Kpp2.ecd' )
plot_file('/home/shafi/wrk/brown/samples/mapk/Kpp_1e-2.ecd' )


#title( figtitle )

#savefig( 'figs/' + figtitle + '.png', dpi=80 )

show()

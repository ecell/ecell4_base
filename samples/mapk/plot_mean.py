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

def resample( x, y, interval ):
    start = x[0]
    end = x[-1]
    tgrid = numpy.mgrid[start:end:interval]
    indices = numpy.searchsorted( x, tgrid )

    return tgrid, y.take( indices )
    

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


    data = load( filename )
    x = data[:,0]
    y = add_columns( data, ycolumns )

    return x, y


def plot_file( filename ):

    x, y = get_data( filename )

    #plot_theory( N_K, N_P, Keq, x[-1] )
    plot( x, y, '-' )

    #psd( y )
    #ylim( 1, 5e4 )


import glob
import os

pattern = sys.argv[1]
globpattern = pattern.replace('ALL','*')

figtitle = os.path.basename( os.path.splitext( pattern )[0] )
print title
#print globpattern
filelist = glob.glob( globpattern )

interval = 30. / 1000

data = []

for filename in filelist:
    x, y = load_data( filename )
    x, y = resample( x, y, interval )
    data.append( y )

y = numpy.array( data )
y = y.mean(0)

plot( x, y )



title( figtitle )

#savefig( 'figs/' + figtitle + '.png', dpi=80 )

show()

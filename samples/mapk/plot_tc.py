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

def add_columns( data, ycolumns ):

    y = numpy.array([ data[:,col] for col in ycolumns ]) 

    y = y.sum(0)

    return y

def get_data( filename ):
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
    #xlim(0,120)


import glob
import os

for filename in sys.argv[1:]:
    plot_file( filename )

#pattern = sys.argv[1]
#globpattern = pattern.replace('ALL','*')

#figtitle = os.path.basename( os.path.splitext( pattern )[0] )
#print title
#print globpattern
#filelist = glob.glob( globpattern )
#print filelist
#for filename in filelist:
#    print filename
#    plot_file( filename )

#title( figtitle )

#savefig( 'figs/' + figtitle + '.png', dpi=80 )

show()

#!/usr/bin/env python


import sys

import numpy
import scipy.io
from matplotlib.pylab import *

from fractionS import *

N_A = 6.0221367e23

E2 = 5
V = 1e-15

def plot_theory( K ):

    N = 1000
    minE1 = 0.1
    maxE1 = 100.
    e1array = numpy.mgrid[minE1:maxE1:(maxE1-minE1)/N]

    farray = [ fraction_Sp( E1, E2, K ) for E1 in e1array ]
    farray = numpy.array( farray )
    #print farray

    semilogx( e1array/E2, farray, label='K = %f' % K )

def file_mean( filename, skip ):
    ycolumns = [2,]
    #ycolumns = [2,6]
    #ycolumns = [3,5]
    #ycolumns = [2,6,3,5]

    data = load( filename )
    x = data[:,0]
    y = data[:,ycolumns[0]]

    start = x.searchsorted( skip )
    if len(x)<=start:
        return None

    x = x[start:]
    y = y[start:]
    #print x[-1]

    xdiff = x[1:] - x[:-1] 
    yscaled = y[:-1] * xdiff
    yscaledmean = yscaled.sum() / ( x[-1] - x[0] )
    print yscaledmean, y.mean()
    #return y.mean()
    return yscaledmean



import glob
import os

S_tot = 300.0

model = 'pushpull'
Keq_str = '0.05'
#Keq_str = '5'
koff_ratio_str = '0.1'
#koff_ratio_str = '0.5'
#koff_ratio_str = '0.9'
#koff_ratio_str = '0'
N_P = 10
V = '1e-14'
T = '400'
#mode = 'normal'
#mode = 'localized'
mode = 'single'

skip = float(T) *0.9

dir = sys.argv[1]
outdir = sys.argv[2]
#pattern = sys.argv[2]
#globpattern = pattern.replace('ALL','*') + '_*.dat'


for N_K in range( 40 ):
    globpattern = \
        string.join( ( model, Keq_str, koff_ratio_str, str(N_K), 
                       str(N_P), V, mode, '*' ), '_' )
    print globpattern
    filelist = glob.glob( dir + os.sep + globpattern )
    if not filelist:
        continue

    data = []

    for file in filelist:
        print file
        res = file_mean( file, skip )
        if res:
            data.append( res )
    data = numpy.array( data )
    print data
    data /= S_tot
    mean = data.mean()
    std_err = data.std()/math.sqrt(len(data))
    print mean, std_err

    errorbar( float(N_K)/N_P, mean, yerr=std_err )

plot_theory( float( Keq_str ) )

figtitle = string.join( ( model, Keq_str, koff_ratio_str, 'ALL', 
                          str(N_P), V, mode ), 
                        '_' )
title( figtitle )

show()
#savefig( outdir + '/' + figtitle + '.png', dpi=80 )


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

    x = x[start:]
    y = y[start:]

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
#Keq_str = '0.05'
Keq_str = '5'
koff_ratio_str = '0.5'
#koff_ratio_str = '0'
N_P = 5
V = '1e-14'
T = '100'
mode = 'normal'
#mode = 'localized'

skip = float(T) / 2

dir = sys.argv[1]
#pattern = sys.argv[2]
#globpattern = pattern.replace('ALL','*') + '_*.dat'


for N_K in ( 1, 2, 3, 4, 5 ):
    globpattern = \
        string.join( ( model, Keq_str, koff_ratio_str, str(N_K), 
                       str(N_P), V, mode, T, '*' ), '_' )
    print globpattern
    filelist = glob.glob( dir + os.sep + globpattern )
    data = numpy.zeros(len(filelist))

    for i, file in enumerate( filelist ):
        print file
        data[i] = file_mean( file, skip )
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

#show()
savefig( 'figs2/' + figtitle + '.png', dpi=80 )


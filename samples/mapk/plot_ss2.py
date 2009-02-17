#!/usr/bin/env python


import sys
import os
import string

import numpy
import scipy.io
from matplotlib.pylab import *

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
    ycolumns = [1,]
    #ycolumns = [2,6]
    #ycolumns = [3,5]
    #ycolumns = [2,6,3,5]

    f = open( filename )
    f.seek( -1000, os.SEEK_END )
    lines = f.readlines()

    lastline = lines[-1]

    lastlinedata = lastline.split()
    if lastlinedata[0] < skip-1:
            raise 'oops'

    y = float( lastlinedata[1] )

    return y

    
#     data = load( filename )
#     x = data[:,0]
#     y = data[:,ycolumns[0]]

#     start = x.searchsorted( skip ) - 1
#     if len(x)<=start:
#         return None

#     return y[start]

#     x = x[start:]
#     y = y[start:]
#     #print x[-1]

#     xdiff = x[1:] - x[:-1] 
#     yscaled = y[:-1] * xdiff
#     yscaledmean = yscaled.sum() / ( x[-1] - x[0] )
#     print yscaledmean, y.mean()
#     #return y.mean()
#     return yscaledmean



import glob
import fnmatch
import os

model = 'mapk4'
V_str = '1e-15'
D_ratio_str = '1'
mode = 'fixed'
N_K_total_str = '300'
#ti_str = '1e-2'
#ti_str = '0'

T = '300'


skip = float(T) #*0.95

dir = sys.argv[1]
#outdir = sys.argv[2]
#pattern = sys.argv[2]
#globpattern = pattern.replace('ALL','*') + '_*.dat'


#os.chdir(dir)

x_all = []
mean_all = []
std_err_all = []


for Kpp_ratio_str in ['0','.3','.5','.7','1']:

    x = []
    mean = []
    std_err = []

    for ti_str in ['0','1e-6','1e-5','1e-4','1e-3','1e-2','1e-1']:

        globpattern = \
            '_'.join( ( model, V_str, D_ratio_str, mode, N_K_total_str,
                        Kpp_ratio_str, ti_str, 'normal',
                            '*' ) ) +\
                            '_tc.dat'

        filelist = glob.glob( dir + os.sep + globpattern )

        if not filelist:
            continue
        #print globpattern

        data = []

        for file in filelist:
            print file
            res = file_mean( file, skip )

            data.append( res )

        data = numpy.array( data )
        data /= int( N_K_total_str )
            
        x.append( float( ti_str ) )
        mean.append( data.mean() )
        std_err.append( data.std()/math.sqrt(len(data)) )

        print x, mean, std_err

    x_all.append( x )
    mean_all.append( mean )
    std_err_all.append( std_err )



axes([.15,.13,.1,.8])
#plot( [1e-6,1], [0,1] )

for i in range( len( x_all ) ):
    errorbar( numpy.array(x_all[i])+1e-18, mean_all[i], yerr=std_err_all[i], 
              fmt='s' )

xlim( [-1e-7,1e-7] )
ylim( [-0.02, 1.01] )

xticks( [0, ], ['0',], size=22 )
yticks( size=22 )

ylabel('[Kpp] / [K]_total', size=28 )

#xscale( 'symlog' )


axes([.26,.13,.7,.8])

#semilogx( [5e-7,1], [0,1] )


for i in range( len( x_all ) ):
    errorbar( numpy.array(x_all[i])+1e-18, mean_all[i], yerr=std_err_all[i], 
              fmt='s' )

xscale( 'log' )


xlim( [1e-7,0.5] )
ylim( [-0.02, 1.01] )


xticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1],['1 us', '10', '100', '1 ms', '10', '100'],size=22)
#xticks( [1e-6, 1e-3, 1e0], ['1 us', '1 ms', '1 s'], size=22 )
yticks( [],[] )


xlabel('t_half', size=28 )

    


show()
#savefig( outdir + '/' + figtitle + '.png', dpi=80 )


#!/usr/bin/env python


import sys
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

model = 'mapk5'
V_str = '1e-15'
D_ratio_str = '1'
#N_KK_str = 
N_P_str = '30'
N_K_total_str = '120'
#ti_str = '1e-2'
ti_str = '0'

T = '60'


skip = float(T) #*0.95

dir = sys.argv[1]
#outdir = sys.argv[2]
#pattern = sys.argv[2]
#globpattern = pattern.replace('ALL','*') + '_*.dat'


#os.chdir(dir)

#for ti_str in ['0','1e-6','1e-4','1e-2']:
for D_ratio_str in ['0.25','0.5','1','2','4']:

    x = []
    mean = []
    std_err = []

    for N_KK in range( 1, 60 ):
        globpattern = \
            '_'.join( ( model, V_str, D_ratio_str, str( N_KK ), '*',
                        N_K_total_str, ti_str,\
                            '*' ) ) +\
                            '_tc.dat'

        filelist = glob.glob( dir + os.sep + globpattern )

        if not filelist:
            continue
        #print globpattern

        for N_P in range( 60 ):

            fnpattern = \
                '_'.join( ( model, V_str, D_ratio_str, str( N_KK ), str( N_P ),
                            N_K_total_str, ti_str,\
                                '*' ) ) +\
                                '_tc.dat'
            filelist2 = fnmatch.filter( filelist, dir + os.sep + fnpattern )
            if not filelist2:
                continue
            #print filelist2


            data = []

            for file in filelist2:
                print file
                res = file_mean( file, skip )

                data.append( res )

            data = numpy.array( data )
            data /= int( N_K_total_str )
            
            x.append( float(N_KK)/float(N_P) )
            mean.append( data.mean() )
            std_err.append( data.std()/math.sqrt(len(data)) )

            print x, mean, std_err

            break


    axes([.135,.13,.8,.8])
    semilogx( x, mean )
    errorbar( x, mean, yerr=std_err, fmt='k+' )

    xlim( [0.005,200] )
    ylim( [-0.02, 1.01] )
    xticks( [1e-2, 1e-1, 1, 1e1, 1e2], ['0.01', '0.1', '1', '10', '100'], size=22 )
    yticks( size=22 )

    xlabel('[KK] / [P]', size=28 )
    ylabel('[Kpp] / [K]_total', size=28 )

    axes([.62,.19,.28,.30])
    plot( x, mean )
    xticks( [.2,.6,1],size=18 )
    yticks( [0,.1,.2,.3,.4],size=18 )
    xlim( [0.001,1] )
    ylim( [-0.02,0.4] )
    #errorbar( x2, mean2, yerr=std_err, fmt='k+' )


show()
#savefig( outdir + '/' + figtitle + '.png', dpi=80 )


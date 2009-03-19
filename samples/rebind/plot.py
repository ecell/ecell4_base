#!/usr/bin/env/python

# varying kf
# python plot.py 05/data/rebind_1_0.1_ALL_t.dat 05/data/rebind_1_1_ALL_t.dat 05/data/rebind_1_10_ALL_t.dat 

# 0.01 didn't run correctly?
# 05/data/rebind_1_0.01_ALL_t.dat 


# varying D
# python plot.py 07/data/rebind_0.1_10_0_ALL_t.dat 07/data/rebind_1_10_0_ALL_t.dat 07/data/rebind_10_10_0_ALL_t.dat


import sys

import numpy
import scipy.io
from matplotlib.pylab import *

#import _gfrd


N_A = 6.0221367e23

N = 1000

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

    n, bins = numpy.histogram(numpy.log10(data), bins=N, new=True)
    n = n.astype(numpy.floating)
    n /= float(len(data))
    x = (10**bins[1:]+ 10**bins[:-1])/2
    dx = (10**bins[1:]- 10**bins[:-1])
    y = n / dx#+1e-10
    print x.shape, y.shape
    print x, y

    print y.sum()
    loglog( x, y )#, label=filename )


if __name__ == '__main__':

    import os
    import glob

    xmin = 1e-9
    xmax = 1e3

    axes([.16,.16,.8,.8])


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


    ka = 0.092e-18
    #1/(1/


    kon = 5.313e-20
    C_B = 16.6e-9

    x = 10 ** numpy.mgrid[-12:3:.1]
    y = kon * C_B * 1000* N_A * numpy.exp( - kon * C_B * 1000 * N_A * x )
    loglog( x, y, 'k--' )

    print x, y


    x = 10 ** numpy.mgrid[-12:-5.7:.1]
    y = 1e3 * x ** (- 1./2.)
    loglog( x, y, 'k:', lw=2 )


    x = 10 ** numpy.mgrid[-5.5:-2:.1]
    y = 1e-3 * x ** (- 3./2.)
    loglog( x, y, 'k:', lw=2 )
    

    text(7e-9,3e7,r'$p \propto \ t^{-1/2}$', size=20)
    text(2e-5,4e4,r'$p \propto \ t^{-3/2}$', size=20)


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

    yticks( [1e-3,1e0, 1e3, 1e6, 1e9], size=20 )

    leg = legend( 
#         # D
         (r'$D=0.1 \ \ {\rm \mu m^2 / s}$',
          r'$D=1 \ \  {\rm \mu m^2 / s}$',
          r'$D=10 \ \  {\rm \mu m^2 / s}$',
          r'Well-stirred ($D=1$)',
        # kf
#         (r'$k_a = 0.017 \ {\rm nM^{-1} s^{-1}}$',
#          r'$k_a = 0.17 \ \ {\rm nM^{-1} s^{-1}}$',
#          r'$k_a = 1.7 \ \ \ \ {\rm nM^{-1} s^{-1}}$',
                   ),
                 loc=1,
                 shadow=True,
                 pad=0.05
                 )
    for l in leg.get_lines():
        l.set_linewidth(1.5)  # the legend line width


    #xlabel( r'$r / \sigma$', fontsize='large' )
    xlabel( r'$t$', size=24 )
    ylabel( r'$p(t)$', size=24 )
    #xlim( 2e-12, 1e2 )
    xlim( 5e-10, 1e2 )
    ylim( 1.1e-6, 2e9 )
    #solline.set_label( r'theory' )
    #legend( handlelen=0.02, pad=0.02,handletextsep=0.01, labelsep=0.001 )
    #grid()
    show()


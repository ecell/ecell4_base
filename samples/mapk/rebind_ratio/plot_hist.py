#!/usr/bin/env python

# D=1

# python plot_hist.py "." mapk3_1e-15_1_fixed_1e-1_normal_ALL_reactions.rebind mapk3_1e-15_1_fixed_1e-2_normal_ALL_reactions.rebind mapk3_1e-15_1_fixed_1e-3_normal_ALL_reactions.rebind mapk3_1e-15_1_fixed_1e-4_normal_ALL_reactions.rebind mapk3_1e-15_1_fixed_1e-5_normal_ALL_reactions.rebind mapk3_1e-15_1_fixed_1e-6_normal_ALL_reactions.rebind mapk3_1e-15_1_fixed_0_normal_ALL_reactions.rebind mapk3_1e-15_1_fixed_0_normal_ALL_reactions.rebind


# t_half = 1e-6

# python plot_hist.py "." mapk3_1e-15_0.25_fixed_1e-6_normal_ALL_reactions.rebind  mapk3_1e-15_0.5_fixed_1e-6_normal_ALL_reactions.rebind  mapk3_1e-15_1_fixed_1e-6_normal_ALL_reactions.rebind  mapk3_1e-15_2_fixed_1e-6_normal_ALL_reactions.rebind  mapk3_1e-15_4_fixed_1e-6_normal_ALL_reactions.rebind 




from matplotlib.pylab import *

import math
import numpy

import sys
import re

def plot_hist( filename, xmin, xmax, N, pattern=None ):

    file = open( filename )

    data=[]

    for line in file.readlines():
        line = line.split()
        t = float(line[0])
        eventType = line[1]
        if t == 0:
            print 'skip zero'
            continue 
        if pattern == None or pattern.match( eventType ):
            data.append( t )

    data = numpy.array(data)

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

    loglog( 10**bins, n+1e-10 )#, label=filename )



if __name__ == '__main__':


    N=50


    pattern = re.compile( sys.argv[1] )
    
    xmin = 1e-12
    xmax = 9
    
    axes([.12,.13,.8,.8])



    for filename in sys.argv[2:]:

        plot_hist( filename, xmin, xmax, N, pattern )


    xlabel( 'Second phosphorylation times', size=22 )

    xticks( [1e-12, 1e-9, 1e-6, 1e-3, 1], 
            [r'${\rm 1 ps}$',
             r'${\rm 1 ns}$',
             r'${\rm 1 \mu s}$',
             r'${\rm 1 ms}$',
             r'${\rm 1 s}$'],
            size=20 )
    
    xlim( xmin, xmax )
    ylim( 4, 5e3 )


    show()

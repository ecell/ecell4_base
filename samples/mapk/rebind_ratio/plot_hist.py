#!/usr/bin/env python

# D=1

# python plot_hist.py "." mapk3_1e-15_1_fixed_1e-1_normal_ALL_reactions.rebind mapk3_1e-15_1_fixed_1e-2_normal_ALL_reactions.rebind mapk3_1e-15_1_fixed_1e-3_normal_ALL_reactions.rebind mapk3_1e-15_1_fixed_1e-4_normal_ALL_reactions.rebind mapk3_1e-15_1_fixed_1e-5_normal_ALL_reactions.rebind mapk3_1e-15_1_fixed_1e-6_normal_ALL_reactions.rebind mapk3_1e-15_1_fixed_0_normal_ALL_reactions.rebind


# t_half = 1e-6

# python plot_hist.py "." mapk3_1e-15_0.03125_fixed_1e-6_normal_ALL_reactions.rebind  mapk3_1e-15_0.0625_fixed_1e-6_normal_ALL_reactions.rebind  mapk3_1e-15_0.25_fixed_1e-6_normal_ALL_reactions.rebind mapk3_1e-15_1_fixed_1e-6_normal_ALL_reactions.rebind  mapk3_1e-15_4_fixed_1e-6_normal_ALL_reactions.rebind 

# t_half = 1e-2

# python plot_hist.py "." mapk3_1e-15_0.03125_fixed_1e-2_normal_ALL_reactions.rebind  mapk3_1e-15_0.0625_fixed_1e-2_normal_ALL_reactions.rebind  mapk3_1e-15_0.25_fixed_1e-2_normal_ALL_reactions.rebind mapk3_1e-15_1_fixed_1e-2_normal_ALL_reactions.rebind  mapk3_1e-15_4_fixed_1e-2_normal_ALL_reactions.rebind 







from matplotlib.pylab import *

import math
import numpy

import sys
import re

def plot_hist( filename, xmin, xmax, N, pattern=None, factor=1.0 ):

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

    logxmin = math.log10(xmin)
    logxmax = math.log10(xmax)
    
    tick=(logxmax-logxmin)/N
    loggrid = numpy.mgrid[logxmin:logxmax:tick]
    grid = numpy.exp(loggrid)
    print len(data)
    print grid, xmin,xmax
    
    n, bins = numpy.histogram(numpy.log10(data), bins=N, new=True)
    print n
    n = n.astype(numpy.floating)
    n /= float(len(data))
    n *= factor

    #x = 10**bins[:-1]
    x = (10**bins[1:] + 10**bins[:-1]) / 2
    dx = (10**bins[1:]- 10**bins[:-1])
    y = n / dx    #  n+1e-10
    print x, y
    print (y*dx).sum()
    return loglog( x, y  )#, label=filename )


def plot_hist2( filename, xmin, xmax, N, pattern=None, factor=1.0 ):

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

    data.sort()

    i = 0
    p = 5
    x = []
    y = []

    ld = len(data)
    while i+p < ld:
        slice = data[i:i+p]
        min, max = slice.min(), slice.max()
        x.append((min + max) / 2)
        y.append(1.0 / (max - min))
        i += p

    y = numpy.array(y,numpy.floating)
    y /= float(len(data))
    y *= factor

    return loglog( x, y  )#, label=filename )



if __name__ == '__main__':


    import numpy

    N=40


    pattern = re.compile( sys.argv[1] )
    
    #xmin = 1e-12
    xmin = 1e-8
    xmax = 50
    
    axes([.16,.16,.8,.8])

    Dlist = [0.03e-12,0.06e-12,0.25e-12,1e-12, 4e-12]

    lines=[]

    for n, filename in enumerate( sys.argv[2:] ):

        D = Dlist[n]

        #         sigma = 5e-9
        #         kD = 4 * numpy.pi * sigma * D
        #         k_a = 9.2e-20#1.6e9 / (1000*6e23)
        #factor = D * ( 1 + (k_a / kD ) )

        factor = 1
        print 'factor', factor
        line = plot_hist( filename, xmin, xmax, N, pattern, factor )
        lines.append(line)


    xlabel( 'Second phosphorylation times', size=26 )
    ylabel( 'Relative frequency', size=26 )
    #ylabel( r'$p(t) \cdot D ( 1 + (k_a / kD))$', size=26 )

    xticks( [1e-12, 1e-9, 1e-6, 1e-3, 1], 
            [r'${\rm 1 ps}$',
             r'${\rm 1 ns}$',
             r'${\rm 1 \mu s}$',
             r'${\rm 1 ms}$',
             r'${\rm 1 s}$'],
            size=24 )
    yticks( size=18 )
    
    xlim( xmin, xmax )
    ylim( 2e-4, 5e5 )

    leg = legend( lines, (r'$D=0.03 \ \ {\rm \mu m^2 / s}$',
                         r'$D=0.06 \ \  {\rm \mu m^2 / s}$',
#                          #              r'$D=0.13 \ \  {\rm \mu m^2 / s}$',
                          r'$D=0.25 \ \  {\rm \mu m^2 / s}$',
                          r'$D=1.0 \ \  {\rm \mu m^2 / s}$',
                          r'$D=4.0 \ \  {\rm \mu m^2 / s}$',
                         ),
                 loc=3,
                 shadow=True,
                 pad=0.05
                 )
    for l in leg.get_lines():
        l.set_linewidth(1.5)  # the legend line width

    show()

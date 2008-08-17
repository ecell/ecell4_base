#!/usr/bin/env python

#   python plot_hist.py "." mapk3_1e-15_4_fixed_1e-1_normal_ALL_reactions.rebind mapk3_1e-15_4_fixed_1e-2_normal_ALL_reactions.rebind mapk3_1e-15_4_fixed_1e-4_normal_ALL_reactions.rebind mapk3_1e-15_4_fixed_1e-6_normal_ALL_reactions.rebind mapk3_1e-15_4_fixed_0_normal_ALL_reactions.rebind
#    python plot_hist.py "." mapk3_1e-15_2_fixed_1e-1_normal_ALL_reactions.rebind mapk3_1e-15_2_fixed_1e-2_normal_ALL_reactions.rebind mapk3_1e-15_2_fixed_1e-4_normal_ALL_reactions.rebind mapk3_1e-15_2_fixed_1e-6_normal_ALL_reactions.rebind mapk3_1e-15_2_fixed_0_normal_ALL_reactions.rebind
#   python plot_hist.py "." mapk3_1e-15_1_fixed_1e-1_normal_ALL_reactions.rebind mapk3_1e-15_1_fixed_1e-2_normal_ALL_reactions.rebind mapk3_1e-15_1_fixed_1e-4_normal_ALL_reactions.rebind mapk3_1e-15_1_fixed_1e-6_normal_ALL_reactions.rebind mapk3_1e-15_1_fixed_0_normal_ALL_reactions.rebind
#   python plot_hist.py "." mapk3_1e-15_0.5_fixed_1e-1_normal_ALL_reactions.rebind mapk3_1e-15_0.5_fixed_1e-2_normal_ALL_reactions.rebind mapk3_1e-15_0.5_fixed_1e-4_normal_ALL_reactions.rebind mapk3_1e-15_0.5_fixed_1e-6_normal_ALL_reactions.rebind mapk3_1e-15_0.5_fixed_0_normal_ALL_reactions.rebind
#   python plot_hist.py "." mapk3_1e-15_0.25_fixed_1e-1_normal_ALL_reactions.rebind mapk3_1e-15_0.25_fixed_1e-2_normal_ALL_reactions.rebind mapk3_1e-15_0.25_fixed_1e-4_normal_ALL_reactions.rebind mapk3_1e-15_0.25_fixed_1e-6_normal_ALL_reactions.rebind mapk3_1e-15_0.25_fixed_0_normal_ALL_reactions.rebind



# t_half = 1e-6

# python plot_hist.py "." mapk3_1e-15_0.25_fixed_1e-6_normal_ALL_reactions.rebind  mapk3_1e-15_0.5_fixed_1e-6_normal_ALL_reactions.rebind  mapk3_1e-15_1_fixed_1e-6_normal_ALL_reactions.rebind  mapk3_1e-15_2_fixed_1e-6_normal_ALL_reactions.rebind  mapk3_1e-15_4_fixed_1e-6_normal_ALL_reactions.rebind 




from matplotlib.pylab import *

import math
import numpy

import sys
import re

N=50


pattern = re.compile( sys.argv[1] )

xmin = 1e-14
xmax = 10


for filename in sys.argv[2:]:

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
    
    n, bins = numpy.histogram(numpy.log(data), bins=N)

    print bins, n

    loglog( numpy.exp(bins), n+1e-10, clip_on=False )#, label=filename )

xlabel( 'second phosphorylation time [s]' )
#legend()
xlim( xmin, xmax )
ylim( 1, 1e4 )


show()

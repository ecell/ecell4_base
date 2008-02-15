#!/usr/bin/env python


import sys

import numpy
import scipy.io
from matplotlib.pylab import *

from fractionS import *

N_A = 6.0221367e23

E2 = 5
kcat = 7.1
V = 1e-15

def plot_theory( E1, E2, K, maxt ):

    frac = 1.0 - fraction_S( E1, E2, K )
    x = [0.0, maxt]
    y = [frac,frac]

    plot( x, y )


def plot_data( data, xcolumn,  ycolumns, St ):

    x = data[:,xcolumn]
    y = numpy.array([ data[:,col] for col in ycolumns ]) 
    y = y[0] + y[1]
    y /= float( St )
    plot( x, y )


K = 0.1
E1 = 3
St = 300
maxt = 1

for E1 in sys.argv[1:]:
    E1 = int( E1 )
    filename = 'pushpull-0_1-%d_timecourse.dat' % E1
    data = load( filename )
    maxt = data[-1,0]
    plot_data( data, 0, (2,6), St )
    plot_theory( E1, E2, 0.1, maxt )




show()

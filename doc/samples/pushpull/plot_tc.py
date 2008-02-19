#!/usr/bin/env python


import sys

import numpy
import scipy.io
from matplotlib.pylab import *

from fractionS import *

N_A = 6.0221367e23

E2 = 5
V = 1e-15

def plot_theory( E1, E2, K, maxt ):

    frac = fraction_S( E1, E2, K )
    x = [0.0, maxt]
    y = [frac,frac]

    plot( x, y )


def plot_data( data, xcolumn,  ycolumns, St ):

    x = data[:,xcolumn]
    y = numpy.array([ data[:,col] for col in ycolumns ]) 
    y = y[0] + y[1]
    y /= float( St )
    plot( x, y )


K = sys.argv[1]
St = 300

for E1 in sys.argv[2:]:
    E1 = int( E1 )
    filename = 'pushpull-%s-0.2-%d.dat' % ( K, E1 )
    data = load( filename )
    maxt = data[-1,0]
    plot_data( data, 0, (3,5), St )
    plot_theory( E1, E2, float( K ), maxt )




show()

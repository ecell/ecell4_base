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
    print farray

    semilogx( e1array/E2, farray, label='K = %f' % K )


plot_theory( 0.01 )
plot_theory( 0.05 )
plot_theory( 0.1 )
plot_theory( 1 )
plot_theory( 10 )
legend()
show()

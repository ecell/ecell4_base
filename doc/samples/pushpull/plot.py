#!/usr/bin/env python


import sys

import numpy
import scipy.io
from matplotlib.pylab import *

N_A = 6.0221367e23

E2 = 5
St = 500e-9
kcat = 0.15e9
V = 1e-15

def fraction_S( E1, E2, St, kcat, K ):

    V1 = kcat * E1 / ( V * N_A )
    V2 = kcat * E2 / ( V * N_A )

    num = E1 - E2 - ( E1 + E2 ) * K + \
        math.sqrt( (E1 - E2) ** 2 + 2 * K * ( E1 - E2 ) ** 2 + \
                       ( ( E1 + E2 ) * K ) ** 2 )
    den = 2 * ( E1 - E2 )

    return num / den




def plot_theory( Km ):
    N = 1000
    minE1 = 0.1
    maxE1 = 100.
    e1array = numpy.mgrid[minE1:maxE1:(maxE1-minE1)/N]

    farray = [ fraction_S( E1, E2, St, kcat, Km ) for E1 in e1array ]
    farray = numpy.array( farray )
    print farray

    semilogx( e1array/E2, farray, label='theory' )



plot_theory( 0.1 )
plot_theory( 1 )
plot_theory( 10 )
show()

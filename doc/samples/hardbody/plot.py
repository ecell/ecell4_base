#!/usr/bin/env/python

import sys

import numpy
import scipy.io
from matplotlib.pylab import *


N_A = 6.0221367e23

def plot_data( N, T ):

    mean = T.mean()
    std_err = T.std()/math.sqrt(len(T))


T=0.001
V=40e-15
#N=100
#C=
D100 =[0.142389059067,0.113878965378,0.115298986435]
D300 =[0.42008805275,0.552675008774,0.504842042923]
D1000 = [2.62250494957,2.66188097,2.69252300262]
D3000 = [15.1535229683,16.8973779678,14.9991471767]

D10000=[150.92720294,
D30000=[
D100000=[

Ns = [100,300,1000,3000,10000,30000,100000]

data1 = [ Ns,
         [D100,D300,D1000,D3000,D10000,D30000,D100000]]




loglog( data1[0] , data1[1], 'o-', label='Vol. = 1e-15 L' )
#loglog( data2[0] , data2[1], 'o-', label='# particles = 600' )
#loglog( data3[0] , data3[1], 'o-', label='Conc. = 1e-6 M' )

xlabel( '# molecules' )
#xlabel( 'Concentration [M]' )
ylabel( 'real time [sec]' )
legend()
show()




#>>> _gfrd.S_irr( .0001 * 1e-8**2/1e-12, 1e-8, 10 * 1e-8 * 1e-12, 1e-12, 1e-8 )
#0.99116163945434221



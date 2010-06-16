#!/usr/bin/env python

#  PYTHONPATH=../.. python plot.py irr.-2.out 0.000000125  irr.-1.out 0.00000125  irr.0.out 0.0000125 irr.1.out 0.000125 irr.2.out 0.00125 irr.3.out 0.0125

# irr.-3.out 0.0000000125 


import sys

import numpy
import scipy.io
from matplotlib.pylab import *


import _gfrd
from p_irr import p_irr

N_A = 6.0221367e23

N = 1000

sigma = 5e-9
r0 = sigma
D_tot = 2e-12
kf = 100 * sigma *  D_tot
#kf = 0

tau = sigma*sigma / D_tot

rmin = sigma

def load_data(filename):
    infile = open(filename)
    data = array([float(x) for x in infile.read().split()], numpy.float)
    infile.close()

    return data
    
def plot_sol(t):

    rmax = 3.1 * math.sqrt(6 * D_tot * t) + rmin

  
    logrmin = math.log(rmin)
    logrmax = math.log(rmax)
    
    tick=(logrmax-logrmin)/N
    loggrid = numpy.mgrid[logrmin:logrmax:tick]
    grid = numpy.exp(loggrid)

    parray = array([p_irr(r, t, r0, kf, D_tot, sigma) for r in grid])

    return loglog(grid / sigma , parray * sigma, 'k-')[0]
    #plot(rarray / sigma , parray, 'k-', label='theory')

def plot_hist(data, T, i):

    bins = 30

    nonreactions = numpy.compress(data >= sigma, data)
    print 'max', max(nonreactions)
    hist, r = numpy.histogram(numpy.log(nonreactions), 
                              bins=bins)
    r = r[:-1]
    histsum = hist.sum()
    S_sim = float(len(nonreactions)) / len(data)
    print 'S_sim', S_sim
    hist = hist.astype(numpy.float)

    r = numpy.concatenate([r, [r[-1] - r[-2]]])
    r = numpy.exp(r)

    xticks = r[1:]-r[:-1]

    hist /= len(data) * xticks

    r = r[:-1] + (xticks * .5)
    #print 'x', x
    #pStyles = ['o', '^', 'v', '<', '>', 's', '+']
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

    loglog(r / sigma, hist * sigma, colors[i] + 'o', 
           label=r'$T = \tau^{%d}$' % round(math.log10(T/tau)))

    

if __name__ == '__main__':

    axes([.15,.15,.8,.8])


    for i in range(len(sys.argv[1:])/2):
        filename = sys.argv[i*2+1]
        T = float(sys.argv[i*2+2])
        print filename,T
        data = load_data(filename)
        plot_hist(data, T, i)
        solline = plot_sol(T)



    xlabel(r'$r / \sigma$', size=28)
    ylabel(r'$p_{irr}$', size=28)
    xlim(0.9, 2.2e2)
    ylim(2e-6, 2e1)
    xticks([1, 10, 100], ['1', '10', '100'], size=22)
    yticks(size=18)
    #solline.set_label(r'theory')
    #legend(handlelen=0.02, pad=0.02,handletextsep=0.01, labelsep=0.001)
    #grid()
    show()




#>>> _gfrd.S_irr(.0001 * 1e-8**2/1e-12, 1e-8, 10 * 1e-8 * 1e-12, 1e-12, 1e-8)
#0.99116163945434221



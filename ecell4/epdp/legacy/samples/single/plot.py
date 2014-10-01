#!/usr/bin/env/python

'''

'''

import sys

import numpy
import scipy.io
from matplotlib.pylab import *


import gfrdbase

N_A = gfrdbase.N_A

D = 1e-12


def plot_sol(t, rmax):
    rmin = 0
    
    N = 100
    rtick = (rmax - rmin) / N
    rarray = numpy.mgrid[rmin:rmax:rtick]
    
    parray = array([gfrdbase.p_free(r, t, D) for r in rarray])
    
    print rarray, parray

    plot(rarray, parray, 'b-')



def plot_file(infilename, t, maxr):

    bins = 20

    infile = open(infilename)
    data = array([float(x) for x in infile.read().split()], numpy.float)
    infile.close()

    hist, lower_edges = numpy.histogram(data, bins=bins,
                                        range=[0, maxr])
    print 'hist', hist
    
    hist = hist.astype(numpy.float)
    
    xtick = lower_edges[2]-lower_edges[1]
    hist /= len(data) * xtick
    x = lower_edges + (xtick * .5)

    plot(x[: -1], hist, '.', label=infilename)


if __name__ == '__main__':


    for i in range(len(sys.argv[1:])/2):
        filename = sys.argv[i*2+1]
        t = float(sys.argv[i*2+2])

        rmax = 3 * math.sqrt(6 * D * t)
        plot_sol(t, rmax)
        plot_file(filename, t, rmax)
    
    legend()
    show()



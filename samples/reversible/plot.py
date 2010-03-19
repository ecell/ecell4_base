#!/usr/bin/env/python

# PYTHONPATH=../.. python plot.py rev.-3.out p_rev.-3.tsv 0.0000000125 rev.-2.out p_rev.-2.tsv 0.000000125 rev.-1.out p_rev.-1.tsv 0.00000125 rev.0.out p_rev.0.tsv 0.0000125 rev.1.out p_rev.1.tsv 0.000125 rev.2.out p_rev.2.tsv 0.00125 rev.3.out p_rev.3.tsv 0.0125

import sys

import numpy
import scipy.io
from matplotlib.pylab import *

#import _gfrd

infilename = sys.argv[1]


N_A = 6.0221367e23

sigma = 5e-9

#r0 = sigma
D_tot = 2e-12
#kf = 10 * sigma * D

tau = sigma*sigma / D_tot
rmin = sigma


def load_data(filename):
    infile = open(filename)
    data = array([float(x) for x in infile.read().split()], numpy.float)
    infile.close()

    return data
    
def plot_sol(filename, t):

    rmax = 3.1 * math.sqrt(6 * D_tot * t) + rmin

    data = scipy.io.read_array(filename)
    rarray, parray = numpy.transpose(data)
    mask = numpy.less_equal(rarray, rmax)
    rarray = numpy.compress(mask, rarray)
    parray = numpy.compress(mask, parray)

    return loglog(rarray / sigma, parray * sigma, 'k-')[0]



def plot_hist(data, T, i):

    bins = 30

    nonreactions = numpy.compress(data >= sigma, data)
    print 'max', max(nonreactions)
    hist, r = numpy.histogram(numpy.log(nonreactions), 
                              bins=bins)
    r = r[:-1]  # new numpy.histogram returns len(r)=len(hist)+1
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

    for i in range(len(sys.argv[1:])/3):
        simfilename = sys.argv[i*3+1]
        solfilename = sys.argv[i*3+2]
        T = float(sys.argv[i*3+3])
        print simfilename,solfilename,T
        data = load_data(simfilename)
        plot_hist(data, T, i)
        solline = plot_sol(solfilename, T)


    xlabel(r'$r / \sigma$', size=28)
    ylabel(r'$p_{rev}$', size=28)

    xlim(0.9, 2.2e2)
    ylim(2e-6, 2e1)
    xticks([1, 10, 100], ['1', '10', '100'], size=22)
    yticks(size=18)
    solline.set_label(r'theory')
    #legend(handlelen=0.02, pad=0.02,handletextsep=0.01, labelsep=0.001)
    #grid()
    show()


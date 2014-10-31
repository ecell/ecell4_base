#!/usr/bin/env python


import sys
import string
import fnmatch

import numpy
import scipy.io
from matplotlib.pylab import *



from fractionS import *

N_A = 6.0221367e23

E2 = 5

def plot_theory(K):

    N = 1000
    minE1 = 0.1
    maxE1 = 100.
    e1array = numpy.mgrid[minE1:maxE1:(maxE1-minE1)/N]

    farray = [fraction_Sp(E1, E2, K) for E1 in e1array]
    farray = numpy.array(farray)
    #print farray

    semilogx(e1array/E2, farray, label='K = %f' % K)

def file_mean(filename, skip):

    f = open(filename)
    f.seek(-1000, os.SEEK_END)
    lines = f.readlines()

    lastline = lines[-1]

    lastlinedata = lastline.split()
    if lastlinedata[0] < skip-1:
            raise 'oops'
    Sp = int(lastlinedata[5])
    PSp = int(lastlinedata[6])
    print Sp, PSp
    y = float(Sp + PSp)
    print lastlinedata
    return y

# def file_mean(filename, skip):
#     ycolumns = [2, ]
#     #ycolumns = [2,6]
#     #ycolumns = [3,5]
#     #ycolumns = [2,6,3,5]

#     data = loadtxt(filename)
#     x = data[:,0]
#     y = data[:,ycolumns[0]]

#     start = x.searchsorted(skip)
#     if len(x)<=start:
#         return None

#     x = x[start:]
#     y = y[start:]
#     #print x[-1]

#     xdiff = x[1:] - x[:-1] 
#     yscaled = y[:-1] * xdiff
#     yscaledmean = yscaled.sum() / (x[-1] - x[0])
#     print yscaledmean, y.mean()
#     #return y.mean()
#     return yscaledmean



import glob
import os

S_tot = 200.0

E_tot = 20

model = 'pushpull'
Keq_str = '0.03'
#Keq_str = '5'
#koff_ratio_str = '0.1'
koff_ratio_str = '0.5'
#koff_ratio_str = '0.1'
#koff_ratio_str = '0'
V = '1e-16'
T = '100'
#mode = 'normal'
#mode = 'localized'
#mode = 'immobile'

skip = float(T) *0.9

dir = sys.argv[1]
outdir = sys.argv[2]
#pattern = sys.argv[2]
#globpattern = pattern.replace('ALL','*') + '_*.dat'


def plot(Keq_str, mode, V):

    N_P=None
    for N_K in range(20):
        globpattern = '_'.join((model, Keq_str, koff_ratio_str, '200', str(N_K), 
                                '*', '*', mode, '*')) + '_tc.dat'
        print globpattern
        filelist = glob.glob(dir + os.sep + globpattern)
        if not filelist:
            continue
            
        N_P = E_tot - N_K

        fnpattern = \
            '_'.join((model, Keq_str, koff_ratio_str, '200', str(N_K), 
                      str(N_P), V, mode, '*')) + '_tc.dat'

        filelist2 = fnmatch.filter(filelist, dir + os.sep + fnpattern)
        if not filelist2:
            continue

        data = []


        for file in filelist2:
            print file
            res = file_mean(file, skip)
            if res:
                data.append(res)
        data = numpy.array(data)
        print data
        data /= S_tot
        mean = data.mean()
        std_err = data.std()/math.sqrt(len(data))
        print mean, std_err

        errorbar(float(N_K)/N_P, mean, yerr=std_err, fmt='k+')

        #if N_K != N_P:
        errorbar(float(N_P)/N_K, 1.0 - mean, yerr=std_err, fmt='k+')

    plot_theory(float(Keq_str))

    xlim(0.02,50)
    ylim(0,1.0)

    figtitle = string.join((model, Keq_str, koff_ratio_str, 'ALL', 
                            V, mode), 
                           '_')
    #title(figtitle)

    #show()
    savefig(outdir + '/' + figtitle + '.png', dpi=80)
    cla()


for Keq in ['0.03', '0.1', '0.3', '1', '3']:
    for V in ['1e-16', '1e-15']:
        for mode in ['localized', 'immobile', 'normal']:
            
            plot(Keq, mode, V)


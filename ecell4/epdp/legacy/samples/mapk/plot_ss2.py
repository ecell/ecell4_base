#!/usr/bin/env python


import sys
import os
import string

import numpy
import scipy.io
from matplotlib.pylab import *

N_A = 6.0221367e23

E2 = 5
V = 1e-15

def load_theory():

    data = load('ss2_ode.dat')

    ti = data[0:len(data):2][:,0]
    data0 = data[0:len(data):2][:,1]
    data1 = data[1:len(data):2][:,1]

    return ti, data0, data1


def file_mean(filename, skip):
    ycolumns = [1, ]
    #ycolumns = [2,6]
    #ycolumns = [3,5]
    #ycolumns = [2,6,3,5]

    f = open(filename)
    f.seek(-1000, os.SEEK_END)
    lines = f.readlines()

    lastline = lines[-1]

    lastlinedata = lastline.split()
    if lastlinedata[0] < skip-1:
            raise 'oops'

    y = float(lastlinedata[1])

    return y

    
#     data = load(filename)
#     x = data[:,0]
#     y = data[:,ycolumns[0]]

#     start = x.searchsorted(skip) - 1
#     if len(x)<=start:
#         return None

#     return y[start]

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
import fnmatch
import os

model = 'mapk4'
V_str = '1e-15'
D_ratio_str = '1'
mode = 'fixed'
N_K_total_str = '300'
#ti_str = '1e-2'
#ti_str = '0'

T = '300'


skip = float(T) #*0.95

#dir = sys.argv[1]
dir = '11/data'
#outdir = sys.argv[2]
#pattern = sys.argv[2]
#globpattern = pattern.replace('ALL','*') + '_*.dat'


#os.chdir(dir)

x_all = []
mean_all = []
std_err_all = []


for Kpp_ratio_str in ['0','.3','.7','1']:

    x = []
    mean = []
    std_err = []

    for ti_str in ['0','1e-6','1e-5','1e-4','1e-3','1e-2','1e-1']:

        globpattern = \
            '_'.join((model, V_str, D_ratio_str, mode, N_K_total_str,
                      Kpp_ratio_str, ti_str, 'normal',
                          '*')) +\
                            '_tc.dat'

        filelist = glob.glob(dir + os.sep + globpattern)

        if not filelist:
            continue
        #print globpattern

        data = []

        for file in filelist:
            print file
            res = file_mean(file, skip)

            data.append(res)

        data = numpy.array(data)
        data /= int(N_K_total_str)
            
        x.append(float(ti_str))
        mean.append(data.mean())
        std_err.append(data.std()/math.sqrt(len(data)))

        print x, mean, std_err

    x_all.append(x)
    mean_all.append(mean)
    std_err_all.append(std_err)


ti, theory0, theory1 = load_theory()

axes([.15,.13,.1,.8])
#plot([1e-6,1], [0,1])

for i in range(len(x_all)):
    errorbar(numpy.array(x_all[i])+1e-18, mean_all[i], yerr=std_err_all[i], 
             fmt='s')

plot(ti[:2],theory0[:2],'k--')
plot(ti[:2],theory1[:2],'k--')

xlim([-1e-7,1e-7])
ylim([-0.02, 1.01])

xticks([0, ], ['$0$', ], size=22)
yticks([0,0.2,0.4,0.6,0.8,1], 
       ['$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1.0$'], size=22)

ylabel(r'$\rm{[Kpp] / [K]_{total}}$', size=28)

#xscale('symlog')


axes([.26,.13,.7,.8])

#semilogx([5e-7,1], [0,1])


for i in range(len(x_all)):
    errorbar(numpy.array(x_all[i])+1e-18, mean_all[i], yerr=std_err_all[i], 
             fmt='s')

semilogx(ti,theory0,'k--')
semilogx(ti,theory1,'k--')

xscale('log')


xlim([1e-7,0.5])
ylim([-0.02, 1.01])


xticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1],
       [r'$1 \mu s$', '$10$', '$100$', r'$1 ms$', '$10$', '$100$'],size=22)
#xticks([1e-6, 1e-3, 1e0], ['1 us', '1 ms', '1 s'], size=22)
yticks([],[])


xlabel(r'${\tau}_{\rm rel}$', size=28)

    

show()
#savefig(outdir + '/' + figtitle + '.png', dpi=80)


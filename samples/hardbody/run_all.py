#!/usr/bin/env python
import sys
import numpy
import math

import run_single

outfile = open(sys.argv[1], 'w')

T = 11.696

REPEAT = 3


def run_set(name, V_list, N_list, T_list):
    
    outfile.write('%s = [\n' % name)
    for i in range(len(V_list)):
        outfile.write('# T=%g, N=%g, V=%g\n' % 
                      (T_list[i], N_list[i], V_list[i]))
        run_times = []
        est_times = []
        for c in range(REPEAT):
            run_time = run_single.run_single(T_list[i], V_list[i], N_list[i])
            est_time = run_time * (T / T_list[i])
            run_times.append(run_time)
            est_times.append(est_time)
        outfile.write('# run_times = %s\n' % str(run_times))
        outfile.write('%s,\n' % str(est_times))
        outfile.flush()
    outfile.write(']\n')




Vv = [40e-15, ] * 10
Nv = [30,100,300,1000,3000,10000,30000,100000,300000,1000000]
#Tv = [1e-0, 1e-1, 1e-1, 1e-2, 1e-3, 1e-3, 1e-4, 1e-4, 1e-5,1e-6]

Tv = [max(1e-6,
          min(T,
              1e4 / math.pow(N, 5.0 / 3.0))) for N in Nv]

Vc = [40e-17, 13e-16, 40e-16, 13e-15, 40e-15, 13e-14, 40e-14, 13e-13, 40e-13,
      13e-12]
Nc = [30,   100,   300, 1000, 3000,10000,30000,100000,300000,1000000]
#Tc = [1e-1, 1e-2, 1e-2, 1e-3, 1e-3, 1e-3, 1e-3, 1e-4, 1e-4, 1e-4]#,1e-5

Tc = [max(1e-3,
          min(T,
              1e1 / math.pow(float(N), 1.0))) for N in Nc]
#Tc = [1e-3] * 10


V300 = [40e-14, 40e-15, 40e-16, 40e-17, 40e-18, 40e-19, 40e-20]#, 40e-21]#,40e-22]
N300 = [300, ] * 7
#T300 = [1e-2, 1e-2, 1e-3, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]#, 1e-8]

T300 = [1e8 / math.pow(1.0/V, 2.0 / 3.0) for V in V300]


V3000 = [40e-13, 40e-14, 40e-15, 40e-16, 40e-17, 40e-18, 40e-19]#, 40e-20]#,40e-21]
N3000 = [3000, ] * 8
#T3000 = [1e-2, 1e-2, 1e-3, 1e-3, 1e-3, 1e-4, 1e-5, 1e-6]#, 1e-7]
T3000 = [1e7 / math.pow(1.0/V, 2.0 / 3.0) for V in V3000]


run_set('data_V', Vv, Nv, Tv); outfile.write('\n\n')
run_set('data_C', Vc, Nc, Tc); outfile.write('\n\n')
run_set('data_N300', V300, N300, T300); outfile.write('\n\n')
run_set('data_N3000', V3000, N3000, T3000); outfile.write('\n\n')


#!/usr/bin/env python
import sys
import numpy
import math

import run_single

prefix = 'out_'

T = 10.

REPEAT = 3


def run_set(outfile, name, V_list, N_list, T_list):
    
    outfile.write('%s = [\n' % name)
    for i in range(len(V_list)):
        outfile.write('# T=%g, N=%g, V=%g\n' % 
                      (T_list[i], N_list[i], V_list[i]))
        run_times = []
        est_times = []
        for c in range(REPEAT):
            run_time, steps, stepspersec = run_single.run_single(T_list[i], 
                                                                 V_list[i], 
                                                                 N_list[i])
            est_time = run_time * (T / T_list[i])
            run_times.append(run_time)
            est_times.append(est_time)
        outfile.write('# steps= %d, steps/sec= %f, steps/N= %f\n'\
                          % (steps, stepspersec, float(steps) / N_list[i]))
        outfile.write('# run_times = %s\n' % str(run_times))
        outfile.write('%s,\n' % str(est_times))
        outfile.flush()
    outfile.write(']\n')


def run_set_bd(outfile, name, V_list, N_list, T_list, dt_factor):
    
    outfile.write('# dt_factor = %g\n' % dt_factor)
    outfile.write('%s = [\n' % name)
    for i in range(len(V_list)):
        outfile.write('# T=%g, N=%g, V=%g\n' % 
                      (T_list[i], N_list[i], V_list[i]))
        run_times = []
        est_times = []
        for c in range(REPEAT):
            run_time, steps, stepspersec = run_single.run_single_bd(T_list[i], 
                                                                    V_list[i], 
                                                                    N_list[i],
                                                                    dt_factor)
            est_time = run_time * (T / T_list[i])
            run_times.append(run_time)
            est_times.append(est_time)
        outfile.write('# steps= %d, steps/sec= %f, steps/N= %f\n'\
                          % (steps, stepspersec, float(steps) / N_list[i]))
        outfile.write('# run_times = %s\n' % str(run_times))
        outfile.write('%s,\n' % str(est_times))
        outfile.flush()
    outfile.write(']\n')




#Vv = [40e-15, ] * 12
Vv = [1e-12, ] * 11
Nv = [100,300,1000,3000,10000,30000,100000,300000,1000000,3000000,10000000]#,30000000]
#Tv = [1e-0, 1e-1, 1e-1, 1e-2, 1e-3, 1e-3, 1e-4, 1e-4, 1e-5,1e-6]

# Tv = [max(1e-3,
#           min(T,
#               3e4 / math.pow(N, 5.0 / 3.0))) for N in Nv]

Tv = [max(1e-5,
          min(T,
              1e0 / math.pow(N, 2.0 / 3.0))) for N in Nv]

#Tv = [1e-3,] * 11

# Vc = [40e-17, 13e-16, 40e-16, 13e-15, 40e-15, 13e-14, 40e-14, 13e-13, 40e-13,
#       13e-12,40e-12,13e-11,40e-12]
Vc = [3.33e-15,1e-14, 3.33e-14,1e-13, 3.33e-13,1e-12, 3.33e-12,1e-11, 3.33e-11,1e-10, 3.33e-10,1e-9]#,3.33e-9]

Nc = [100,      300,   1000,    3000,  10000,  30000,100000,300000,1000000,3000000,10000000,30000000]#,100000000]
#Tc = [1e-1, 1e-2, 1e-2, 1e-3, 1e-3, 1e-3, 1e-3, 1e-4, 1e-4, 1e-4]#,1e-5

Tc = [max(1e-3,
          min(T,
              1e1 / math.pow(float(N), 1.0))) for N in Nc]
#Tc = [1e-3] * 10


V300 = [40e-14, 40e-15, 40e-16, 40e-17, 40e-18, 40e-19, 40e-20, 13.3e-20]#,40e-22]
N300 = [300, ] * 8
#T300 = [1e-2, 1e-2, 1e-3, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]#, 1e-8]
T300 = [1e7 / math.pow(1.0/V, 2.0 / 3.0) for V in V300]


V3000 = [40e-13, 40e-14, 40e-15, 40e-16, 40e-17, 40e-18, 40e-19, 13.3e-19]#,40e-21]
N3000 = [3000, ] * 8
#T3000 = [1e-2, 1e-2, 1e-3, 1e-3, 1e-3, 1e-4, 1e-5, 1e-6]#, 1e-7]
T3000 = [2e6 / math.pow(1.0/V, 2.0 / 3.0) for V in V3000]


VBD = [1e-12] * 7
NBD = [100,300,1000,3000,10000,30000,100000]
#NBD = [100,1000,100000]
TBD = [1e-4 / N for N in NBD]
#TBD = [1e-5 / N for N in NBD]

TBD2 = [1e-3 / N for N in NBD]


VBD300 = V300
NBD300 = N300
TBD300 = [5e-7] * len(VBD300)

BD_DTFACTOR = 1e-5

if __name__ == '__main__':
    mode = sys.argv[1]
    outfile = open(prefix+mode+'.py','w'); 
    dataname = 'data_' + mode
    if mode == 'V':
        run_set(outfile, dataname, Vv, Nv, Tv); outfile.write('\n\n')
    elif mode == 'C':
        run_set(outfile, dataname, Vc, Nc, Tc); outfile.write('\n\n')
    elif mode == 'N300':
        run_set(outfile, dataname, V300, N300, T300); outfile.write('\n\n')
    elif mode == 'N3000':
        run_set(outfile, dataname, V3000, N3000, T3000); outfile.write('\n\n')
    elif mode == 'BD':
        run_set_bd(outfile, dataname, VBD, NBD, TBD, BD_DTFACTOR); outfile.write('\n\n')

    elif mode == 'BD300':
        run_set_bd(outfile, dataname, VBD300, NBD300, TBD300, BD_DTFACTOR); outfile.write('\n\n')
    # elif mode == 'BD2':
    #     run_set_bd(outfile, dataname, VBD, NBD, TBD2, 1e-4); outfile.write('\n\n')

    # just for large # particles stress tests
    elif mode == 'NE6':
        run_set(outfile, dataname, [3.3e-9], [1e6], [1e-3]); outfile.write('\n\n')
    elif mode == 'N3E6':
        run_set(outfile, dataname, [1e-10], [3e6], [1e-3]); outfile.write('\n\n')
    elif mode == 'NE7':
        run_set(outfile, dataname, [3.3e-10], [1e7], [1e-3]); outfile.write('\n\n')
    elif mode == 'NE5BD':
        run_set_bd(outfile, dataname, [1e-12], [1e5], [1e-9], 1e-5); outfile.write('\n\n')



    else:
        raise 'invalid argument'


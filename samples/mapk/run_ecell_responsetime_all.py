import os
import math
import numpy

MODEL_FILE0 = 'model4-0.eml'
MODEL_FILE1 = 'model4.eml'

ESS_FILE = 'run_ecell_tc.py'

DURATION = 1000

# Register jobs.

jobs = {}



N_P = 30
N_KK = 30

N_KPP = 0
N_K = 120
        
cwd = os.getcwd()

for ti_str in ['1e-2', '1e-6']:

    ti = float(ti_str)

    if ti == 0:
        MODEL_FILE = MODEL_FILE0
        KI='invalid'
    else:
        MODEL_FILE = MODEL_FILE1
        KI = math.log(2) / ti


    for D_str in ['%.3g'% 10 ** e for e in numpy.mgrid[-2:2.1:.1]]:

        D = float(D_str)

        OUTFILE = '%s/Kpp_ODE_%s_%s.ecd' % (cwd,D_str,ti_str)

        parameterDict = { 'MODEL_FILE': MODEL_FILE,
                          'DURATION': DURATION,
                          'OUTFILE': OUTFILE,
                          'D': D_str,
                          'N_KK': N_KK, 'N_P': N_P,
                          'N_KPP': N_KPP, 'N_K': N_K, 'KI': KI }
        
        jobID = registerEcellSession(ESS_FILE, parameterDict, [MODEL_FILE, ])
        #jobs[jobID] = [float(N_KK)/N_P, kpp_ratio]

run()

#import sys

# for jobID in jobs.keys():

#     #print " --- job id = %s ---" % jobID
#     sys.stdout.write( '[%s, %s],' % ( jobs[jobID][0],
#                                   getStdout(jobID) ) )

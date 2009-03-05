import os
import math
import numpy

MODEL_FILE0 = 'model4-0.eml'
MODEL_FILE1 = 'model4.eml'

MODEL_FILE1 = 'model4-processive.eml'

ESS_FILE = 'run_ecell_model4.py'

N_KK_P_total = 60
N_K_total = 120

# Register jobs.

jobs = {}


N_KK_list = range(60)

ti_str = '1e-2'
#ti_str = '1e-6'
#ti_str = '0'

ti = float( ti_str )

if ti == 0:
    MODEL_FILE = MODEL_FILE0
    KI='invalid'
else:
    MODEL_FILE = MODEL_FILE1
    KI = math.log( 2 ) / ti

kpp_ratio = 0.5
for N_KK in N_KK_list:

    N_P = N_KK_P_total - N_KK

    N_KPP = N_K_total * kpp_ratio
    N_K = N_K_total - N_KPP

    parameterDict = { 'MODEL_FILE': MODEL_FILE,
                      'N_KK': N_KK, 'N_P': N_P,
                      'N_KPP': N_KPP, 'N_K': N_K, 'KI': KI }
        
    jobID = registerEcellSession( ESS_FILE, parameterDict, [ MODEL_FILE, ])
    jobs[ jobID ] = [ float(N_KK)/N_P, kpp_ratio ]

run()

import sys

for jobID in jobs.keys():

    #print " --- job id = %s ---" % jobID
    sys.stdout.write( '[%s, %s],' % ( jobs[jobID][0],
                                  getStdout( jobID ) ) )

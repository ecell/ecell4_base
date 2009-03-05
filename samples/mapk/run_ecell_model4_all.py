import os
import math
import numpy

MODEL_FILE0 = 'model4-0.eml'
MODEL_FILE1 = 'model4.eml'
ESS_FILE = 'run_ecell_model4.py'

N_K_total = 300

# Register jobs.

jobs = {}



for ti_str in [ '0', '1e-7', '3e-7',
                '1e-6', '3e-6',
                '1e-5', '3e-5',
                '1e-4', '3e-4',
                '1e-3', '3e-3',
                '1e-2', '3e-2',
                '1e-1', '3e-1',
                '1e-0' ]:
    ti = float( ti_str )

    if ti == 0:
        MODEL_FILE = MODEL_FILE0
    else:
        MODEL_FILE = MODEL_FILE1
        KI = math.log( 2 ) / ti


    for kpp_ratio in [ 0.11, 0.66 ]:
        N_KPP = N_K_total * kpp_ratio
        N_K = N_K_total - N_KPP

        if ti != 0:
            parameterDict = { 'MODEL_FILE': MODEL_FILE,
                              'N_KPP': N_KPP, 'N_K': N_K, 'KI': KI }
        else:
            parameterDict = { 'MODEL_FILE': MODEL_FILE,
                              'N_KPP': N_KPP, 'N_K': N_K }
        
        jobID = registerEcellSession( ESS_FILE, parameterDict, [ MODEL_FILE, ])
        jobs[ jobID ] = [ ti_str, kpp_ratio ]

run()

import sys

for jobID in jobs.keys():

    #print " --- job id = %s ---" % jobID
    sys.stdout.write( '%s %s' % ( jobs[jobID][0],
                                  getStdout( jobID ) ) )

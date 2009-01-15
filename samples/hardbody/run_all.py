#!/usr/bin/env python
import sys
import numpy

import run_single

outfile = open( sys.argv[1], 'w' )

T = 1.1696


def run_set( name, V_list, N_list, T_list ):
    
    outfile.write( '%s = [\n' % name )
    for i in range( len( V_list ) ):
        outfile.write( '# T=%g, N=%g, V=%g\n' % 
                       ( T_list[i], N_list[i], V_list[i] ) )
        outfile.write( '[' )
        for c in range( 3 ):
            runtime = run_single.run_single( T_list[i], V_list[i], N_list[i] )
            runtime *= T/T_list[i]
            outfile.write( '%g,' % runtime )
        outfile.write( '],\n' )
        outfile.flush()
    outfile.write( ']\n' )




Vv = [40e-15,] * 10
Nv = [30,100,300,1000,3000,10000,30000,100000,300000,1000000]
Tv = [1e-0, 1e-1, 1e-1, 1e-2, 1e-3, 1e-3, 1e-4, 1e-4, 1e-5,1e-6]




Vc = [40e-17, 13e-16, 40e-16, 13e-15, 40e-15, 13e-14, 40e-14, 13e-13, 40e-13,
      13e-12]
Nc = [30,   100,   300, 1000, 3000,10000,30000,100000,300000,1000000]
Tc = [1e-1, 1e-2, 1e-2, 1e-3, 1e-3, 1e-3, 1e-3, 1e-4, 1e-4, 1e-4]#,1e-5



V300 = [40e-14, 40e-15, 40e-16, 40e-17, 40e-18, 40e-19, 40e-20, 10e-21,40e-21]
N300 = [300,] * 9
T300 = [1e-2, 1e-2, 1e-3, 1e-3, 1e-4, 1e-5, 1e-6, 1e-6, 1e-7]


V3000 = [40e-13, 40e-14, 40e-15, 40e-16, 40e-17, 40e-18, 40e-19, 10e-20,40e-21]
N3000 = [3000,] * 9
T3000 = [1e-2, 1e-2, 1e-3, 1e-3, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]


run_set( 'data_V', Vv, Nv, Tv ); outfile.write( '\n\n' )
run_set( 'data_C', Vc, Nc, Tc ); outfile.write( '\n\n' )
run_set( 'data_N300', V300, N300, T300 ); outfile.write( '\n\n' )
run_set( 'data_N3000', V3000, N3000, T3000 ); outfile.write( '\n\n' )


#!/usr/bin/env python

from egfrd import *
from bd import *

from logger import *
import sys
import time


def run_single( T, V, N ):

    print 'T =', T, '; V= ', V, '; N=', N
    

    L = math.pow( V * 1e-3, 1.0 / 3.0 )

    s = EGFRDSimulator()
    #s = BDSimulator()
    s.setWorldSize( L )

    matrixSize = min( max( 3, int( (3 * N) ** (1.0/3.0) ) ), 120 )
    print 'matrixSize=', matrixSize
    s.setMatrixSize( matrixSize )
    
    box1 = CuboidalSurface( [0,0,0],[L,L,L] )

    D = 1e-12

    A = Species( 'A', D, 2.5e-9 )
    s.addSpecies( A )
    
    s.throwInParticles( A, N, box1 )
    print 'stir'

    stirTime = T * .1
    while 1:
        s.step()
        nextTime = s.getNextTime()
        if nextTime > stirTime:
            s.stop( stirTime )
            break
    print 'reset'
    s.reset()
    print 'reset finish'

    print 'run'
    start = time.time()
    while s.t < T:
        s.step()
        #l.log()

    end = time.time()
    print 'TIMING:\n', end - start

    return end - start


if __name__ == '__main__':
    
    T = float( sys.argv[1] )
    V = float( sys.argv[2] )
    N = int( sys.argv[3] )

    run_single( T, V, N )


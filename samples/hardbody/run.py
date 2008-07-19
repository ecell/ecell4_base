#!/usr/bin/env python

from egfrd import *
from bd import *

from logger import *
import sys
import time

# run
T = float( sys.argv[1] )
V = float( sys.argv[2] )
N = int( sys.argv[3] )

L = math.pow( V * 1e-3, 1.0 / 3.0 )


s = EGFRDSimulator()
#s = BDSimulator()
s.setWorldSize( L )

matrixSize = min( max( 3, int( (3 * N) ** (1.0/3.0) ) ), 100 )
print 'matrixSize=', matrixSize
s.setMatrixSize( matrixSize )

#print int( N ** (1.0/3.0) )
#sys.exit(0)
#s.setMatrixSize(10)
#s.setMatrixSize( 30 )

box1 = CuboidalSurface( [0,0,0],[L,L,L] )

A = Species( 'A', 1e-12, 2.5e-9 )
s.addSpecies( A )

# no reaction

s.throwInParticles( A, N, box1 )
print 'stir'
#stirTime = 0
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
#l = Logger( s, 'hardbody' )
#l.setParticleOutput( ('A',) )
#l.setInterval( 1e-1 )
#l.log()

print 'run'
start = time.time()
while s.t < T:
    s.step()
    #l.log()

end = time.time()
print 'T =', T, '; V= ', V, '; N=', N
print 'TIMING:\n', end - start

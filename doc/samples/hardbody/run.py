#!/usr/bin/env python

from egfrd import *

from logger import *
import sys
import time

# run
T = float( sys.argv[1] )
V = float( sys.argv[2] )
N = int( sys.argv[3] )

L = math.pow( V * 1e-3, 1.0 / 3.0 )

s = EGFRDSimulator( 'normal' )
s.setWorldSize( L )

s.setMatrixSize( max( 3, int( N ** (1.0/3.0) ) ) )
#s.setMatrixSize( 30 )

box1 = CuboidalSurface( [0,0,0],[L,L,L] )

A = Species( 'A', 1e-12, 2.5e-9 )
s.addSpecies( A )

# no reaction

s.throwInParticles( A, N, box1 )
print 'stir'
stirTime = 1e-7
while 1:
    s.step()
    nextTime = s.scheduler.getTopTime()
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

#!/usr/bin/env python

from egfrd import *

from logger import *
import sys
import time

# run
T = 0.01
V = 40e-15

# 1-1.  C=1e-7M
#N=60
# 2.56697297096
# 1-2.  C=2e-7M
#N=120
# 9.34459114075
# 1-3.  C=1e-6M
#N=600
# 296.172708988
# 1-4.  C=2e-6M
T=0.01
N=1000

# 2
#N=600
# 2-1.  C=1e-8M
#V=1e-13 
#5.03202104568
#5.85098314285
#6.20103216171
#5.03884100914

# 2-2.  C=1e-7M
# V=1e-14
#  33.6928880215
#  39.1034328938

# 2-3.  C=1e-6M
# V=1e-15

# 2-4.  C=1e-5M
#V=1e-16
#T=0.001
#446.12724494
#308.7816581739

#3
#C=1e-6M
#3.1
#V=1e-17
#N=6
# 1.37475013733
# 1.96632885933
# 2.04573392868
# 1.54014611244

#3.2
#V=1e-16
#N=60
# 18.4902999401
# 29.3690078259
# 28.4070279598
# 17.4782500267

#3.3
#V=1e-15
#N=600
#279.476380825

#3.4
#V=1e-14
#N=6000
#T=0.001
#1326.545048
#1248.54476881

#3.5
# V=1e-13
# N=60000

# T=0.0001


L = math.pow( V * 1e-3, 1.0 / 3.0 )

s = EGFRDSimulator( 'normal' )
s.setWorldSize( L )

s.setMatrixSize( 20 )

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

print 'TIMING:\t', end - start

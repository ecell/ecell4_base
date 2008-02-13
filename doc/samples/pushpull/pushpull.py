#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

def k_a( k ):
    Dpisigma4 = 4 * numpy.pi * Dtot * sigma
    kon = k / 1e3 / N_A
    k_smol = Dpisigma4
    print 'kon ', kon, 'k_smol', k_smol
    if kon > k_smol:
        print 'ERROR: kon > k_smol, implying k_a is negative'
        sys.exit( 1 )
    ka = 1 / ( ( 1 / kon ) - ( 1 / k_smol ) )
    return ka

def k_d( koff, kon ):
    return k_a( kon ) * koff / kon * N_A * 1e3

def C2N( c, V ):
    return round( c * V * N_A )


radius = 5e-9
sigma = radius * 2
D = 1.0e-12

#V = 2000 * sigma * sigma * sigma * 1000 * 500 #liter
V = 1e-15
L = math.pow( V * 1e-3, 1.0 / 3.0 )

s = EGFRDSimulator()
s.setCellSize( L )
print V, L

print C2N( 498e-9, V )

#sys.exit(0)


box1 = CuboidalSurface( [0,0,0],[L,L,L] )
# not supported yet
#s.addSurface( box1 )


S = Species( 'S', D, radius )
s.addSpecies( S )
P = Species( 'P', D, radius )
s.addSpecies( P )
K = Species( 'K', D, radius )
s.addSpecies( K )
KS = Species( 'KS', D, radius )
s.addSpecies( KS )
Sp = Species( 'Sp', D, radius )
s.addSpecies( Sp )
PSp = Species( 'PSp', D, radius )
s.addSpecies( PSp )

Dtot = D + D

#  1 2 S + K  <-> KS
#  3   KS      -> K + Sp
#  4 5 Sp + P <-> PSp
#  6   PSp     -> P + S

k1 = k2 = 6.0
kb = 1.5
kf = 0.15e9

k_a_f = k_a( kf )
print k_a_f
print k_d( kb, kf )

#sys.exit(0)

s.throwInParticles( P, 5, box1 )
s.throwInParticles( K, 5, box1 )
s.throwInParticles( Sp, 150, box1 )
s.throwInParticles( S, 150, box1 )



# Stir before actually start the sim.

stirTime = 1e-3
while 1:
    s.step()
    nextTime = s.scheduler.getNextTime()
    if nextTime > stirTime:
        s.stop( stirTime )
        break

s.reset()


r1 = BindingReactionType( S, K, KS, k_a_f )
s.addReactionType( r1 )
r2 = UnbindingReactionType( KS, S, K, k_d( kb, kf ) )
s.addReactionType( r2 )
r3 = UnbindingReactionType( KS, K, Sp, k_d( k1, kf ) )
s.addReactionType( r3 )
r4 = BindingReactionType( Sp, P, PSp, k_a_f )
s.addReactionType( r4 )
r5 = UnbindingReactionType( PSp, Sp, P, k_d( kb, kf ) )
s.addReactionType( r5 )
r6 = UnbindingReactionType( PSp, P, S, k_d( k2, kf ) )
s.addReactionType( r6 )



l = Logger( s, 'pushpull-0_05' )
#l.setParticleOutput( ('Ea','X','EaX','Xp','Xpp','EaI') )
l.setInterval( 1e-3 )
l.log()


while s.t < 10:
    s.step()
    #s.dumpPopulation()
    l.log()
    


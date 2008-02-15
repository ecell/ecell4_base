#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

from fractionS import *

def k_a( k ):
    Dpisigma4 = 4 * numpy.pi * Dtot * sigma
    kon = k / 1e3 / N_A
    k_D = Dpisigma4
    print 'kon ', kon, 'k_D', k_D
    if kon > k_D:
        print 'ERROR: kon > k_D.'
        sys.exit( 1 )
    ka = 1 / ( ( 1 / kon ) - ( 1 / k_D ) )
    return ka

def k_d( koff, kon ):
    return k_a( kon ) * koff / kon * N_A * 1e3

def C2N( c, V ):
    return round( c * V * N_A )

Keq = float( sys.argv[1] )
koff_ratio = float( sys.argv[2] )
N_K = int( sys.argv[3] )


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

N_P = 5

fracS = fraction_S( N_K, N_P, Keq )

S_tot = 300
S_conc = S_tot / N_A / V   # in M

N_S = S_tot * fracS
N_Sp = S_tot - N_S

kon = 0.15e9

Keq_S = Keq * S_conc

kcatkoff = Keq_S * kon

koff = kcatkoff * koff_ratio
kcat = kcatkoff - koff


kf = k_a( kon )

print 'kon', kon, 'koff', koff, 'kcat', kcat

#sys.exit(0)

s.throwInParticles( K, N_K, box1 )
s.throwInParticles( P, N_P, box1 )
s.throwInParticles( Sp, N_Sp, box1 )
s.throwInParticles( S, N_S, box1 )

# Stir before actually start the sim.

stirTime = 1e-8
while 1:
    s.step()
    nextTime = s.scheduler.getNextTime()
    if nextTime > stirTime:
        s.stop( stirTime )
        break

s.reset()

#  1 2 S + K  <-> KS
#  3   KS      -> K + Sp
#  4 5 Sp + P <-> PSp
#  6   PSp     -> P + S


r1 = BindingReactionType( S, K, KS, kf )
s.addReactionType( r1 )
r2 = UnbindingReactionType( KS, S, K, k_d( koff, kon ) )
s.addReactionType( r2 )
r3 = UnbindingReactionType( KS, K, Sp, kcat )
s.addReactionType( r3 )
r4 = BindingReactionType( Sp, P, PSp, kf )
s.addReactionType( r4 )
r5 = UnbindingReactionType( PSp, Sp, P, k_d( koff, kon ) )
s.addReactionType( r5 )
r6 = UnbindingReactionType( PSp, P, S, kcat )
s.addReactionType( r6 )



l = Logger( s, 'pushpull-%s-%s-%s' % ( sys.argv[1], sys.argv[2], sys.argv[3] ),
            comment = '# kon=%f; koff=%f; kcat=%f; V=%f; S_tot=%f; N_P=%f' %
            ( kon, koff, kcat, V, S_tot, N_P ) )
#l.setParticleOutput( ('Ea','X','EaX','Xp','Xpp','EaI') )
#l.setInterval( 1e-3 )
l.log()


while s.t < 30:
    s.step()
    #s.dumpPopulation()
    l.log()
    


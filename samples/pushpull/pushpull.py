#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

from fractionS import *


# Args:
# Keq
# koff_ratio
# N_K
# N_P
# V (liter)
# mode:  'normal' 'immobile' 'localized' 'single' 'clustered'
# T

Keq_str = sys.argv[1]
koff_ratio_str = sys.argv[2]
N_S_total = int( sys.argv[3] )
N_K = int( sys.argv[4] )
N_P = int( sys.argv[5] )
V_str = sys.argv[6]
mode = sys.argv[7]
T_str = sys.argv[8]


Keq = float( Keq_str )
koff_ratio = float( koff_ratio_str )
V = float( V_str )
T = float( T_str )

radius = 2.5e-9
sigma = radius * 2
D1 = 1.0e-12



if mode == 'normal':
    D2 = D1
elif mode == 'immobile' or mode == 'localized' or mode == 'single':
    D2 = 0
else:
    raise 'invalid mode'


L = ( V * 1e-3 ) ** ( 1.0 / 3.0 )


s = EGFRDSimulator()
s.setWorldSize( L )

N = N_S_total * 1.1
matrixSize = min( max( 3, int( (3 * N) ** (1.0/3.0) ) ), 60 )
print 'matrixSize=', matrixSize
s.setMatrixSize( matrixSize )

#s.setDtFactor( 1e-5 )

print V, L

print C2N( 498e-9, V )

#sys.exit(0)


box1 = CuboidalSurface( [0,0,0],[L,L,L] )
plain1 = CuboidalSurface( [0,0,0],[0,L,L] )
plain2 = CuboidalSurface( [L/2,0,0],[L/2,L,L] )
# not supported yet
#s.addSurface( box1 )


S = Species( 'S', D1, radius )
s.addSpecies( S )
P = Species( 'P', D2, radius )
s.addSpecies( P )
K = Species( 'K', D2, radius )
s.addSpecies( K )
KS = Species( 'KS', D2, radius )
s.addSpecies( KS )
Sp = Species( 'Sp', D1, radius )
s.addSpecies( Sp )
PSp = Species( 'PSp', D2, radius )
s.addSpecies( PSp )

#fracS = fraction_S( N_K, N_P, Keq )
fracS = 1


S_conc = N_S_total / V * 1e3   # in #/m^3

N_S = N_S_total * fracS
N_Sp = N_S_total - N_S

Dtot = D1 + D2

#ka = k_a( kon, k_D( Dtot, sigma ) )
#ka = 9e9 / N_A / 1e3 # 1/M s -> m^3/s
kD = k_D( Dtot, sigma )

kon = Mtom3( 0.02e9 )
ka = k_a( kon, kD )

Keq_S = Keq * S_conc

kcatkoff = Keq_S * kon
koff = kcatkoff * koff_ratio
kcat = kcatkoff - koff

if mode == 'single':
    kcat1 = kcat * float( N_K ) / float( N_P )
    koff1 = kcatkoff - kcat1
    kcat2 = kcat
    koff2 = koff
else:
    kcat1 = kcat2 = kcat
    koff1 = koff2 = koff


kd1 = k_d( koff, kon, kD )
kd2 = k_d( koff2, kon, kD )

print 'ka', ka, 'kD', kD, 'kd1', kd1, 'kd2', kd2
print 'kon m^3/s', kon, '1/M s', kon * N_A * 1e3
print 'koff1 1/s ', koff1
print 'kcat1 1/s ', kcat1
print 'koff2 1/s ', koff2
print 'kcat2 1/s ', kcat2

assert koff2 >= 0

#sys.exit(0)

if mode == 'normal' or mode == 'immobile':
    s.throwInParticles( K, N_K, box1 )
    s.throwInParticles( P, N_P, box1 )
elif mode == 'localized':
    s.throwInParticles( K, N_K, plain1 )
    s.throwInParticles( P, N_P, plain2 )
elif mode == 'single':
    x = L/2
    yz = L/2
    tl = L/4
    s.placeParticle( K, [ tl, tl, tl ] )
    s.placeParticle( K, [ tl, tl, yz+tl ] )
    s.placeParticle( K, [ tl, yz+tl, tl ] )
    s.placeParticle( K, [ tl, yz+tl, yz+tl ] )
    s.placeParticle( P, [ x+tl, tl, tl ] )
    s.placeParticle( P, [ x+tl, tl, yz+tl ] )
    s.placeParticle( P, [ x+tl, yz+tl, tl ] )
    s.placeParticle( P, [ x+tl, yz+tl, yz+tl ] )
else:
    assert False



s.throwInParticles( Sp, N_Sp, box1 )
s.throwInParticles( S, N_S, box1 )

# Stir before actually start the sim.

stirTime = 1e-7
while 1:
    s.step()
    nextTime = s.scheduler.getTopTime()
    if nextTime > stirTime:
        s.stop( stirTime )
        break

s.reset()

#  1 2 S + K  <-> KS
#  3   KS      -> K + Sp
#  4 5 Sp + P <-> PSp
#  6   PSp     -> P + S


r1 = BindingReactionType( S, K, KS, ka )
s.addReactionType( r1 )
r2 = UnbindingReactionType( KS, S, K, kd1 )
s.addReactionType( r2 )
r3 = UnbindingReactionType( KS, K, Sp, kcat1 )
s.addReactionType( r3 )
r4 = BindingReactionType( Sp, P, PSp, ka )
s.addReactionType( r4 )
r5 = UnbindingReactionType( PSp, Sp, P, kd2 )
s.addReactionType( r5 )
r6 = UnbindingReactionType( PSp, P, S, kcat2 )
s.addReactionType( r6 )



model = 'pushpull'

# 'pushpull-Keq-koff_ratio-N_K-N_P-V-mode.dat'
l = Logger( s, 
            logname = model + '_' + '_'.join( sys.argv[1:8] ) + '_' +\
                sys.env[ 'SGE_TASK_ID' ],
            comment = '@ model=\'%s\'; Keq=%s; koff_ratio=%s\n' %
            ( model, Keq_str, koff_ratio_str ) +
            '#@ V=%s; N_K=%s; N_P=%s; mode=\'%s\'; T=%s\n' % 
            ( V_str, N_K, N_P, mode, T_str ) +
            '#@ kon=%g; koff1=%g; koff2=%g; N_S_total=%s\n' %
            ( kon, koff1, koff2, N_S_total ) +
            '#@ kcat1=%g; kcat2=%g\n' %
            ( kcat1, kcat2 ) +
            '#@ ka=%g; kd1=%g; kd2=%g\n' %
            ( ka, kd1, kd2 ) )
#l.setParticleOutput( ('K','P') )
#l.setInterval( 1e-3 )
#l.writeParticles()
l.log()


while s.t < T:
    s.step()

    if s.lastReaction:
        #log.info( s.dumpPopulation() )
        l.log()
    


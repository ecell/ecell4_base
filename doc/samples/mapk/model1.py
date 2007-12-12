#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

import math

V = 1e-15 #liter
#L = 1e-6
#V = 3.33e-17
L = math.pow( V * 1e-3, 1.0 / 3.0 )

s = EGFRDSimulator()
s.setCellSize( L )

box1 = CuboidalSurface( [0,0,0],[L,L,L] )
# not supported yet
#s.addSurface( box1 )

#D = 2e-12 # run1
#D = 1e-12 # run2
D = 5e-13 # run3
#D = 0.25e-12 # run4

K = Species( 'K', D, 5e-9 )
s.addSpecies( K )
KK = Species( 'KK', D, 5e-9 )
s.addSpecies( KK )
P = Species( 'P', D, 5e-9 )
s.addSpecies( P )
Kp = Species( 'Kp', D, 5e-9 )
s.addSpecies( Kp )
Kpp = Species( 'Kpp', D, 5e-9 )
s.addSpecies( Kpp )
K_KK = Species( 'K_KK', D, 5e-9 )
s.addSpecies( K_KK )
Kp_KK = Species( 'Kp_KK', D, 5e-9 )
s.addSpecies( Kp_KK )
Kpp_KK = Species( 'Kpp_KK', D, 5e-9 )
s.addSpecies( Kpp_KK )
Kpp_P = Species( 'Kpp_P', D, 5e-9 )
s.addSpecies( Kpp_P )
Kp_P = Species( 'Kp_P', D, 5e-9 )
s.addSpecies( Kp_P )

#  1 2   K + KK   <-> K_KK
#  3     K_KK       -> Kp + KK
#  4 5   Kp + KK  <-> Kp_KK
#  6     Kp_KK      -> Kpp + KK 
#  7 8   Kpp + P <-> Kpp_P
#  9     Kpp_P     -> Kp + P
# 10 11  Kp + P  <-> Kp_P
# 12     Kp_P      -> K + P
# 12     Kp_P      -> K + P
# 13     Kpp     -> Kp
# 14     Kp      -> K


Dtot = D + D
sigma = 5e-9 * 2
Dpisigma4 = 4 * numpy.pi * Dtot * sigma

def k_a( k ):
    kon = k / 1e3 / N_A
    k_smol = Dpisigma4
    print kon, k_smol
    return 1 / ( ( 1 / kon ) - ( 1 / k_smol ) )

def k_d( koff, kon ):
    return k_a( kon ) * koff / kon * N_A * 1e3
    

print k_a( 0.02e9)
print k_d( 1.0, 0.02e9 )
#raise ''

def C2N( c ):
    return round( c * V * N_A )

s.throwInParticles( K, C2N( 200e-9 ), box1 )

s.throwInParticles( KK, C2N( 50e-9 ), box1 )
s.throwInParticles( P, C2N( 50e-9 ), box1 )

print ''
while s.t < 1:
    s.step()
s.reset()

r1 = BindingReactionType( K, KK, K_KK, k_a(0.02e9) )
s.addReactionType( r1 )
r2 = UnbindingReactionType( K_KK, K, KK, k_d( 1.0, 0.02e9 ) )
s.addReactionType( r2 )
r3 = UnbindingReactionType( K_KK, Kp, KK, k_d( 1.5, 0.02e9 ) )
s.addReactionType( r3 )

r4 = BindingReactionType( Kp, KK, Kp_KK, k_a( 0.032e9 ) )
s.addReactionType( r4 )
r5 = UnbindingReactionType( Kp_KK, Kp, KK, k_d( 1.0, 0.032e9 ) )
s.addReactionType( r5 )
r6 = UnbindingReactionType( Kp_KK, Kpp, KK, k_d( 15.0, 0.032e9 ) )
s.addReactionType( r6 )

r7 = BindingReactionType( Kpp, P, Kpp_P, k_a( 0.02e9 ) )
s.addReactionType( r7 )
r8 = UnbindingReactionType( Kpp_P, Kpp, P, k_d( 1.0, 0.02e9 ) )
s.addReactionType( r8 )
r9 = UnbindingReactionType( Kpp_P, Kp, P, k_d( 1.5, 0.02e9 ) )
s.addReactionType( r9 )

r10 = BindingReactionType( Kp, P, Kp_P, k_a( 0.032e9 ))
s.addReactionType( r10 )
r11 = UnbindingReactionType( Kp_P, Kp, P, k_d( 1.0, 0.032e9 ) )
s.addReactionType( r11 )
r12 = UnbindingReactionType( Kp_P, K, P, k_d( 15.0, 0.032e9 ) )
s.addReactionType( r12 )

#r13 = UnimolecularReactionType( Kpp, Kp, 1e-1 )
#s.addReactionType( r13 )
#r14 = UnimolecularReactionType( Kp, K, 1e-1 )
#s.addReactionType( r14 )


l = Logger( s, 'mapk3' )
#l.setParticleOutput( ('Ea','X','EaX','Xp','Xpp','EaI') )
l.setInterval( 1e-0 )
l.log()

while s.t < 30:
    s.step()
    s.dumpPopulation()
    l.log()


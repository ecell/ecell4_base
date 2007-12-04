#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

s = EGFRDSimulator()
s.setCellSize( 1e-6 )

box1 = CuboidalSurface( [0,0,0],[1e-6,1e-6,1e-6] )
# not supported yet
#s.addSurface( box1 )

K = Species( 'K', 2e-11, 5e-8 )
s.addSpecies( K )
KK = Species( 'KK', 2e-11, 5e-8 )
s.addSpecies( KK )
P = Species( 'P', 2e-11, 5e-8 )
s.addSpecies( P )
Kp = Species( 'Kp', 2e-11, 5e-8 )
s.addSpecies( Kp )
Kpp = Species( 'Kpp', 2e-11, 5e-8 )
s.addSpecies( Kpp )
K_KK = Species( 'K_KK', 2e-11, 5e-8 )
s.addSpecies( K_KK )
Kp_KK = Species( 'Kp_KK', 2e-11, 5e-8 )
s.addSpecies( Kp_KK )
Kpp_KK = Species( 'Kpp_KK', 2e-11, 5e-8 )
s.addSpecies( Kpp_KK )
Kpp_P = Species( 'Kpp_P', 2e-11, 5e-8 )
s.addSpecies( Kpp_P )
Kp_P = Species( 'Kp_P', 2e-11, 5e-8 )
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

r1 = BindingReactionType( K, KK, K_KK, 1e9 / N_A )
s.addReactionType( r1 )
r2 = UnbindingReactionType( K_KK, K, KK, 1e3 )
s.addReactionType( r2 )
r3 = UnbindingReactionType( K_KK, Kp, KK, 1e3 )
s.addReactionType( r3 )

r4 = BindingReactionType( Kp, KK, Kp_KK, 1e9 / N_A )
s.addReactionType( r4 )
r5 = UnbindingReactionType( Kp_KK, Kp, KK, 1e3 )
s.addReactionType( r5 )
r6 = UnbindingReactionType( Kpp_KK, Kpp, KK, 1e3 )
s.addReactionType( r6 )

r7 = BindingReactionType( Kpp, P, Kpp_P, 1e9 / N_A )
s.addReactionType( r7 )
r8 = UnbindingReactionType( Kpp_P, Kpp, P, 1e3 )
s.addReactionType( r8 )
r9 = UnbindingReactionType( Kpp_P, Kp, P, 1e3 )
s.addReactionType( r9 )

r10 = BindingReactionType( Kp, P, Kp_P, 1e9 / N_A )
s.addReactionType( r10 )
r11 = UnbindingReactionType( Kp_P, Kp, P, 1e3 )
s.addReactionType( r11 )
r12 = UnbindingReactionType( Kp_P, K, P, 1e3 )
s.addReactionType( r12 )

r13 = UnimolecularReactionType( Kpp, Kp, 1e2 )
s.addReactionType( r13 )
r14 = UnimolecularReactionType( Kp, K, 1e2 )
s.addReactionType( r14 )

s.throwInParticles( K, 60, box1 )
s.throwInParticles( KK, 30, box1 )


l = Logger( s, 'dimer' )
#l.setParticleOutput( ('Ea','X','EaX','Xp','Xpp','EaI') )
l.setInterval( 1e-3 )
l.log()


while s.t < 100:
    s.step()
    s.dumpPopulation()

#    l.log()
    


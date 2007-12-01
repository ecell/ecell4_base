#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

s = EGFRDSimulator()
s.setCellSize( 1e-6 )

box1 = CuboidalSurface( [0,0,0],[1e-5,1e-5,1e-5] )
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
KKK = Species( 'KKK', 2e-11, 5e-8 )
s.addSpecies( KKK )
KpKK = Species( 'KpKK', 2e-11, 5e-8 )
s.addSpecies( KpKK )
KppP = Species( 'KppP', 2e-11, 5e-8 )
s.addSpecies( KppP )
KpP = Species( 'KpP', 2e-11, 5e-8 )
s.addSpecies( KpP )

#  1 2   K + KK   <-> KKK
#  3     KKK       -> Kp + KK
#  4 5   Kp + KK  <-> KpKK
#  6     KpKK      -> Kpp + KK 
#  7 8   Kpp + P <-> KppP
#  9     KppP     -> Kp + P
# 10 11  Kp + P  <-> KpP
# 12     KpP      -> K + P

r1 = BindingReactionType( , X, EaX, 1e9 / N_A )
s.addReactionType( r1 )
r2 = UnbindingReactionType( EaX, Ea, X, 1e3 )
s.addReactionType( r2 )
r3 = UnbindingReactionType( EaX, Ea, Xp, 1e3 )
s.addReactionType( r3 )
r4 = BindingReactionType( Ea, Xp, EaXp, 1e9 / N_A )
s.addReactionType( r4 )
r5 = UnbindingReactionType( EaXp, Ea, Xp, 1e3 )
s.addReactionType( r5 )
r6 = UnbindingReactionType( EaXp, Ea, Xpp, 1e3 )
s.addReactionType( r6 )
r7 = BindingReactionType( Ea, Xpp, EaXpp, 1e9 / N_A )
s.addReactionType( r7 )
r8 = UnbindingReactionType( EaXpp, EaI, Xpp, 1e3 )
s.addReactionType( r8 )
r9 = UnimolecularReactionType( EaI, Ea, 1e2 )
s.addReactionType( r9 )

s.throwInParticles( Ea, 60, box1 )
s.throwInParticles( X, 60, box1 )


l = Logger( s, 'dimer' )
#l.setParticleOutput( ('Ea','X','EaX','Xp','Xpp','EaI') )
l.setInterval( 1e-3 )
l.log()


while s.t < 100:
    s.step()
    s.dumpPopulation()

#    l.log()
    


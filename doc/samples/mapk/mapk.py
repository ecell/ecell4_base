#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

s = EGFRDSimulator()
s.setCellSize( 1e-5 )

box1 = CuboidalSurface( [0,0,0],[1e-5,1e-5,1e-5] )
# not supported yet
#s.addSurface( box1 )

Ea = Species( 'Ea', 2e-11, 5e-8 )
s.addSpecies( Ea )
EaX = Species( 'EaX', 2e-11, 5e-8 )
s.addSpecies( EaX )
EaXp = Species( 'EaXp', 2e-11, 5e-8 )
s.addSpecies( EaXp )
EaXpp = Species( 'EaXpp', 2e-11, 5e-8 )
s.addSpecies( EaXpp )
X = Species( 'X', 2e-11, 5e-8 )
s.addSpecies( X )
Xp = Species( 'Xp', 2e-11, 5e-8 )
s.addSpecies( Xp )
Xpp = Species( 'Xpp', 2e-11, 5e-8 )
s.addSpecies( Xpp )
EaI = Species( 'EaI', 2e-11, 5e-8 )
s.addSpecies( EaI )

#  1 2 Ea + X <--> EaX
#  3   EaX -> Ea + Xp
#  4 5 Ea + Xp <-> EaXp
#  6   EaXp -> Ea + Xpp
#  7   Ea + Xpp -> EaXpp
#  8   EaXpp -> EaI + Xpp
#  9   EaI -> Ea.

r1 = BindingReactionType( Ea, X, EaX, 1e9 / N_A )
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
    


#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

s = EGFRDSimulator()
s.setCellSize( 1e-6 )

box1 = CuboidalSurface( [0,0,0],[1e-5,1e-5,1e-5] )
# not supported yet
#s.addSurface( box1 )

Fus3 = Species( 'Fus3', 2e-11, 5e-8 )
s.addSpecies( Fus3 )
Ste7 = Species( 'Ste7', 2e-11, 5e-8 )
s.addSpecies( Ste7 )
Msg5 = Species( 'Msg5', 2e-11, 5e-8 )
s.addSpecies( Msg5 )
Fus3p = Species( 'Fus3p', 2e-11, 5e-8 )
s.addSpecies( Fus3p )
Fus3pp = Species( 'Fus3pp', 2e-11, 5e-8 )
s.addSpecies( Fus3pp )
Fus3Ste7 = Species( 'Fus3Ste7', 2e-11, 5e-8 )
s.addSpecies( Fus3Ste7 )
Fus3pSte7 = Species( 'Fus3pSte7', 2e-11, 5e-8 )
s.addSpecies( Fus3pSte7 )
Fus3ppMsg5 = Species( 'Fus3ppMsg5', 2e-11, 5e-8 )
s.addSpecies( Fus3ppMsg5 )
Fus3pMsg5 = Species( 'Fus3pMsg5', 2e-11, 5e-8 )
s.addSpecies( Fus3pMsg5 )

#  1 2   Fus3 + Ste7   <-> Fus3Ste7
#  3     Fus3Ste7       -> Fus3p + Ste7
#  4 5   Fus3p + Ste7  <-> Fus3pSte7
#  6     Fus3pSte7      -> Fus3pp + Ste7 
#  7 8   Fus3pp + Msg5 <-> Fus3ppMsg5
#  9     Fus3ppMsg5     -> Fus3p + Msg5
# 10 11  Fus3p + Msg5  <-> Fus3pMsg5
# 12     Fus3pMsg5      -> Fus3 + Msg5

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
    


#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

s = EGFRDSimulator()
s.setCellSize( 1e-6 )

box1 = CuboidalSurface( [0,0,0],[1e-6,1e-6,1e-6] )
# not supported yet
#s.addSurface( box1 )

S = Species( 'S', 2e-11, 5e-8 )
s.addSpecies( S )
P = Species( 'P', 2e-11, 5e-8 )
s.addSpecies( P )
K = Species( 'K', 2e-11, 5e-8 )
s.addSpecies( K )
KS = Species( 'KS', 2e-11, 5e-8 )
s.addSpecies( KS )
Sp = Species( 'Sp', 2e-11, 5e-8 )
s.addSpecies( Sp )
PSp = Species( 'PSp', 2e-11, 5e-8 )
s.addSpecies( PSp )

#  1 2 S + K  <-> KS
#  3   KS      -> K + Sp
#  4 5 Sp + P <-> PSp
#  6   PSp     -> P + S

r1 = BindingReactionType( S, K, KS, 1e9 / N_A )
s.addReactionType( r1 )
r2 = UnbindingReactionType( KS, S, K, 1e3 )
s.addReactionType( r2 )
r3 = UnbindingReactionType( KS, K, Sp, 1e3 )
s.addReactionType( r3 )
r4 = BindingReactionType( Sp, P, PSp, 1e9 / N_A )
s.addReactionType( r4 )
r5 = UnbindingReactionType( PSp, Sp, P, 1e3 )
s.addReactionType( r5 )
r6 = UnbindingReactionType( PSp, P, S, 1e3 )
s.addReactionType( r6 )

s.placeParticle( K, [0,0,0] )
s.placeParticle( P, [5e-7,0,0] )

s.throwInParticles( S, 50, box1 )


l = Logger( s, 'pushpull' )
#l.setParticleOutput( ('Ea','X','EaX','Xp','Xpp','EaI') )
l.setInterval( 1e-3 )
l.log()


while s.t < 100:
    s.step()
    s.dumpPopulation()

#    l.log()
    


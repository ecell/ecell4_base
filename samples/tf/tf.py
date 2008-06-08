#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

s = EGFRDSimulator()
s.setWorldSize( 1e-6 )

box1 = CuboidalSurface( [0,0,0],[1e-6,1e-6,1e-6] )
# not supported yet
#s.addSurface( box1 )

O = Species( 'O', 0, 1e-8 )
s.addSpecies( O )
R = Species( 'R', 1e-12, 1e-8 )
s.addSpecies( R )
P = Species( 'R', 1e-12, 1e-8 )
s.addSpecies( P )
OR = Species( 'OR', 0, 1e-8 )
s.addSpecies( OR )
ORp = Species( 'ORp', 0, 1e-8 )
s.addSpecies( ORp )
ORpa = Species( 'ORpa', 0, 1e-8 )
s.addSpecies( ORpa )
T = Species( 'T', 1e-12, 1e-8 )
s.addSpecies( T )
M = Species( 'M', 1e-12, 1e-8 )
s.addSpecies( M )
Mribo = Species( 'Mribo', 1e-12, 1e-8 )
s.addSpecies( Mribo )
#EMPTY = Species( 'EMPTY', 2e-12, 5e-8 )
#s.addSpecies( EMPTY )

#  1 2 O + R <-> OR
#  3 4 O     <-> ORp
#  5   ORp    -> ORpa
#  6   ORpa   -> T + O
#  7   M      -> EMPTY
#  8   M      -> M + Mribo
#  9   Mribo  -> P
# 10   P      -> EMPTY


k_fR = 6e9 * 1000 / N_A
k_bR = 0.1  # 1 - 0.01
k_fRp = 38
k_bRp = 0.5
k_OC = 1 # 0.3 - 3
t_clear = 1  # should not be poisson
t_elon = 50 # 50-100
k_dm = 0.019
k_ribo = 5 * k_dm
k_dp = 2.4e-4
t_trans = 30


r1 = BindingReactionType( O, R, OR, k_fR )
s.addReactionType( r1 )
r2 = UnbindingReactionType( OR, O, R, k_bR )
s.addReactionType( r2 )
r3 = UnimolecularReactionType( O, ORp, k_fRp )
s.addReactionType( r3 )
r4 = UnimolecularReactionType( ORp, O, k_bRp )
s.addReactionType( r4 )
r5 = UnimolecularReactionType( ORp, ORpa, k_OC )
s.addReactionType( r5 )
r6 = UnbindingReactionType( ORpa, T, O, 1/t_clear )
s.addReactionType( r6 )
r7 = DecayReactionType( M, k_dm )
s.addReactionType( r7 )
r8 = UnbindingReactionType( M, M, Mribo, k_ribo )
s.addReactionType( r8 )
r9 = UnimolecularReactionType( Mribo, P, 1/t_trans )
s.addReactionType( r9 )
r10 = DecayReactionType( P, k_dp )
s.addReactionType( r10 )

s.placeParticle( O, [0,0,0] )

#s.throwInParticles( R, 50, box1 )


l = Logger( s, 'pushpull' )
#l.setParticleOutput( ('Ea','X','EaX','Xp','Xpp','EaI') )
l.setInterval( 1e-3 )
l.log()


while s.t < 1000:
    s.step()
    s.dumpPopulation()

#    l.log()
    


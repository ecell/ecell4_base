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


r1 = BindingReactionRule( O, R, OR, k_fR )
s.addReactionRule( r1 )
r2 = UnbindingReactionRule( OR, O, R, k_bR )
s.addReactionRule( r2 )
r3 = UnimolecularReactionRule( O, ORp, k_fRp )
s.addReactionRule( r3 )
r4 = UnimolecularReactionRule( ORp, O, k_bRp )
s.addReactionRule( r4 )
r5 = UnimolecularReactionRule( ORp, ORpa, k_OC )
s.addReactionRule( r5 )
r6 = UnbindingReactionRule( ORpa, T, O, 1/t_clear )
s.addReactionRule( r6 )
r7 = DecayReactionRule( M, k_dm )
s.addReactionRule( r7 )
r8 = UnbindingReactionRule( M, M, Mribo, k_ribo )
s.addReactionRule( r8 )
r9 = UnimolecularReactionRule( Mribo, P, 1/t_trans )
s.addReactionRule( r9 )
r10 = DecayReactionRule( P, k_dp )
s.addReactionRule( r10 )

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
    


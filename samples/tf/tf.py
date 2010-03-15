#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

w = World(1e-6, 3)
s = EGFRDSimulator(w)

box1 = CuboidalRegion( [0,0,0],[1e-6,1e-6,1e-6] )

m = ParticleModel()
O = m.new_species_type( 'O', 0, 1e-8 )
R = m.new_species_type( 'R', 1e-12, 1e-8 )
P = m.new_species_type( 'R', 1e-12, 1e-8 )
OR = m.new_species_type( 'OR', 0, 1e-8 )
ORp = m.new_species_type( 'ORp', 0, 1e-8 )
ORpa = m.new_species_type( 'ORpa', 0, 1e-8 )
T = m.new_species_type( 'T', 1e-12, 1e-8 )
M = m.new_species_type( 'M', 1e-12, 1e-8 )
Mribo = m.new_species_type( 'Mribo', 1e-12, 1e-8 )
#EMPTY = m.new_species_type( 'EMPTY', 2e-12, 5e-8 )

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


r1 = createBindingReactionRule( O, R, OR, k_fR )
m.network_rules.add_reaction_rule( r1 )
r2 = createUnbindingReactionRule( OR, O, R, k_bR )
m.network_rules.add_reaction_rule( r2 )
r3 = createUnimolecularReactionRule( O, ORp, k_fRp )
m.network_rules.add_reaction_rule( r3 )
r4 = createUnimolecularReactionRule( ORp, O, k_bRp )
m.network_rules.add_reaction_rule( r4 )
r5 = createUnimolecularReactionRule( ORp, ORpa, k_OC )
m.network_rules.add_reaction_rule( r5 )
r6 = createUnbindingReactionRule( ORpa, T, O, 1/t_clear )
m.network_rules.add_reaction_rule( r6 )
r7 = createDecayReactionRule( M, k_dm )
m.network_rules.add_reaction_rule( r7 )
r8 = createUnbindingReactionRule( M, M, Mribo, k_ribo )
m.network_rules.add_reaction_rule( r8 )
r9 = createUnimolecularReactionRule( Mribo, P, 1/t_trans )
m.network_rules.add_reaction_rule( r9 )
r10 = createDecayReactionRule( P, k_dp )
m.network_rules.add_reaction_rule( r10 )

s.setModel( m )

s.placeParticle( O, [0,0,0] )

#s.throwInParticles( R, 50, box1 )


l = Logger( s, 'pushpull' )
#l.setParticleOutput( ('Ea','X','EaX','Xp','Xpp','EaI') )
l.setParticleOutInterval( 1e-3 )
l.log()


while s.t < 1000:
    s.step()
    s.dumpPopulation()

#    l.log()
    


#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

s = EGFRDSimulator()
s.setWorldSize( 1e-5 )

box1 = CuboidalSurface( [0,0,0],[1e-5,1e-5,1e-5] )
# not supported yet
#s.addSurface( box1 )

m = ParticleModel()

Ea = m.new_species_type( 'Ea', 2e-11, 5e-8 )
EaX = m.new_species_type( 'EaX', 2e-11, 5e-8 )
EaXp = m.new_species_type( 'EaXp', 2e-11, 5e-8 )
EaXpp = m.new_species_type( 'EaXpp', 2e-11, 5e-8 )
X = m.new_species_type( 'X', 2e-11, 5e-8 )
Xp = m.new_species_type( 'Xp', 2e-11, 5e-8 )
Xpp = m.new_species_type( 'Xpp', 2e-11, 5e-8 )
EaI = m.new_species_type( 'EaI', 2e-11, 5e-8 )

#  1 2 Ea + X <--> EaX
#  3   EaX -> Ea + Xp
#  4 5 Ea + Xp <-> EaXp
#  6   EaXp -> Ea + Xpp
#  7   Ea + Xpp -> EaXpp
#  8   EaXpp -> EaI + Xpp
#  9   EaI -> Ea.

r1 = createBindingReactionRule( Ea, X, EaX, 1e9 / N_A )
m.network_rules.add_reaction_rule( r1 )
r2 = createUnbindingReactionRule( EaX, Ea, X, 1e3 )
m.network_rules.add_reaction_rule( r2 )
r3 = createUnbindingReactionRule( EaX, Ea, Xp, 1e3 )
m.network_rules.add_reaction_rule( r3 )
r4 = createBindingReactionRule( Ea, Xp, EaXp, 1e9 / N_A )
m.network_rules.add_reaction_rule( r4 )
r5 = createUnbindingReactionRule( EaXp, Ea, Xp, 1e3 )
m.network_rules.add_reaction_rule( r5 )
r6 = createUnbindingReactionRule( EaXp, Ea, Xpp, 1e3 )
m.network_rules.add_reaction_rule( r6 )
r7 = createBindingReactionRule( Ea, Xpp, EaXpp, 1e9 / N_A )
m.network_rules.add_reaction_rule( r7 )
r8 = createUnbindingReactionRule( EaXpp, EaI, Xpp, 1e3 )
m.network_rules.add_reaction_rule( r8 )
r9 = createUnimolecularReactionRule( EaI, Ea, 1e2 )
m.network_rules.add_reaction_rule( r9 )

s.setModel(m)

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
    


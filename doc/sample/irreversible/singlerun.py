#!/usr/bin/env python

from egfrd import *

import sys

s = EGFRDSimulator()
s.setCellSize( 1e-3 )

s.setMaxShellSize( 1e-6 )

A = Species( 'A', 0.0, 5e-8 )
s.addSpecies( A )
B = Species( 'B', 1e-11, 5e-8 )
s.addSpecies( B )
C = Species( 'C', 0.0, 5e-8 )
s.addSpecies( C )

N = Species( 'N', 0.0, 5e-8 )
s.addSpecies( N )

r1 = BindingReactionType( A, B, C, 1e5 / N_A )
s.addReactionType( r1 )

s.setAllRepulsive()

s.placeParticle( A, [0,0,0] )
s.placeParticle( B, [1e-7+1e-18,0,0] )
#s.placeParticle( N, NOWHERE )

while s.t < 100:
    s.step()
    s.dumpPopulation()
    if s.isPopulationChanged:
        print 'reaction'
        sys.exit(1)
    

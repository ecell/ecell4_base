#!/usr/bin/env python

#from egfrd import *
from bd import *

from logger import *
import sys

#s = EGFRDSimulator()

s = BDSimulator()
s.setWorldSize( 1e-6 )
s.setMatrixSize( 10 )

box1 = CuboidalSurface( [0,0,0],[1e-6,1e-6,1e-6] )
# not supported yet
#s.addSurface( box1 )

S = Species( 'S', 1.5e-12, 5e-9 )
s.addSpecies( S )
P = Species( 'P', 1e-12, 7e-9 )
s.addSpecies( P )

r1 = BindingReactionType( S, S, P, 1e7 / N_A )
s.addReactionType( r1 )
r2 = UnbindingReactionType( P, S, S, 1e3 )
s.addReactionType( r2 )

s.throwInParticles( S, 0, box1 )
s.throwInParticles( P, 1000, box1 )

l = Logger( s, 'dimer' )
l.setParticleOutput( ('P','S') )
l.setInterval( 1e-3 )
l.log()


'''
while s.t < 100:
    s.step()
    s.dumpPopulation()
#    l.log()
'''    


def profrun():
    for i in range( 10 ):
        s.step()

import profile
profile.run('profrun()', 'fooprof')
import pstats
pstats.Stats('fooprof').sort_stats('time').print_stats(30)


sys.exit(1)


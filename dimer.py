#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

s = EGFRDSimulator()
s.setCellSize( 1e-5 )


box1 = CuboidalSurface( [0,0,0],[1e-5,1e-5,1e-5] )
# not supported yet
#s.addSurface( box1 )

S = Species( 'S', 2e-11, 5e-8 )
s.addSpecies( S )
P = Species( 'P', 1e-11, 7e-8 )
s.addSpecies( P )

r1 = BindingReactionType( S, S, P, 1e10 / N_A )
s.addReactionType( r1 )
r2 = UnbindingReactionType( P, S, S, 1e3 )
s.addReactionType( r2 )

s.throwInParticles( S, 0, box1 )
s.throwInParticles( P, 60, box1 )

l = Logger( s, 'dimer' )
l.setParticleOutput( ('P','S') )
l.setInterval( 1e-3 )
l.log()


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
'''

#!/usr/bin/env python

from egfrd import *
from bd import *

from logger import *
import sys

s = EGFRDSimulator()
#s = BDSimulator()

L = 5e-6
#L = 2e-6
#L = 5e-8
#L = 8e-7
s.setWorldSize( L )
s.setMatrixSize( 20 )

box1 = CuboidalSurface( [0,0,0], [L,L,L] )
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

s.throwInParticles( S, 1500, box1 )
s.throwInParticles( P, 1500, box1 )

l = Logger( s, 'dimer' )
l.setParticleOutput( ('P','S') )
l.setInterval( 1e-3 )
l.log()



#while s.t < 100:
#    s.step()

#s.dumpPopulation()
#l.log()



def profrun():
    #while s.stepCounter < 6000:
    for _ in range( 12000 ):
        s.step()
        #logging.info( s.dumpPopulation() )


try:
    import cProfile as profile
except:
    import profile
profile.run('profrun()', 'fooprof')

import pstats
pstats.Stats('fooprof').sort_stats('time').print_stats(40)


sys.exit(1)


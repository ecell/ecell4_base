#!/usr/bin/env python

from egfrd import *
from bd import *
from surface import *
from gfrdbase import *
from logger import *
import sys

s = EGFRDSimulator()
#s = BDSimulator()

#L = 5e-6
L = 2e-6
#L = 5e-8
#L = 3e-7
s.setWorldSize( L )
s.setMatrixSize( 10 )

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

s.throwInParticles( S, 150, box1 )
s.throwInParticles( P, 150, box1 )

l = Logger( s, 'dimer' )
l.setParticleOutInterval( 1e-7 )
l.log()


#while s.t < 100:
#    s.step()

#s.dumpPopulation()
#l.log()



def profrun():
    #while s.stepCounter < 6000:
    for _ in xrange( 15000 ):
        s.step()
        #l.log()
        #logging.info( s.dumpPopulation() )


try:
    import cProfile as profile
except:
    import profile
profile.run('profrun()', 'fooprof')
#profrun()

import pstats
pstats.Stats('fooprof').sort_stats('time').print_stats(40)


sys.exit(1)


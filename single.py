#!/usr/bin/env python

from gfrd import *

from logger import *
import sys

s = GFRDSimulator()
s.setCellSize( 1e-5 )
S = Species( 'S', 2e-11, 5e-8 )
s.addSpecies( S )
P = Species( 'P', 0.0, 7e-8 )
s.addSpecies( P )
SP = Species( 'SP', 0.0, 7e-8 )
s.addSpecies( SP )

r1 = BindingReactionType( S, P, SP, 1e18 / N_A )
s.addReactionType( r1 )
#r2 = UnimolecularReactionType( P, S, 1e3 )
r2 = UnbindingReactionType( SP, S, P, 5e1 )
s.addReactionType( r2 )

#s.setAllRepulsive()

s.placeParticle( 'P', [0,0,0] )
s.throwInParticles( 'S', 1 )


l = Logger( s, 'dimer' )
l.setInterval( 1e-6 )
l.log()


while s.t < 5:
    s.step()
    l.log()
    print 't = ', s.t
    

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

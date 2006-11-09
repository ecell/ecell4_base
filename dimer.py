#!/usr/bin/env python

from brown import *

from logger import *
import sys

s = Simulator()
s.setSize( 1e-5 )
S = Species( 'S', 2e-11, 5e-8 )
s.addSpecies( S )
P = Species( 'P', 1e-11, 7e-8 )
s.addSpecies( P )

r1 = BindingReactionType( S, S, P, 1e18 / N_A )
s.addReactionType( r1 )
#r2 = UnimolecularReactionType( P, S, 1e3 )
r2 = UnbindingReactionType( P, S, S, 5e1 )
s.addReactionType( r2 )

s.setAllRepulsive()

s.throwInParticles( 'S', 600 )
s.throwInParticles( 'P', 0 )

l = Logger( s, 'dimer' )
l.interval=1e-6
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

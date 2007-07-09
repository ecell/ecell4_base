#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

s = EGFRDSimulator()
s.setCellSize( 1e-5 )


box1 = CuboidalSurface( [0,0,0],[1e-5,1e-5,1e-5] )
# not supported yet
#s.addSurface( box1 )

P = Species( 'P', 2e-11, 1e-7 )
s.addSpecies( P )

s.setAllRepulsive()

s.throwInParticles( 'P', 600, box1 )

l = Logger( s, 'simple' )
l.setParticleOutput( ('P',) )
l.setInterval( 1e-3 )
l.log()

while s.t < 1:
    s.step()
    l.log()
    

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

#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

s = EGFRDSimulator()
size = 5e-6
s.setCellSize( size )


box1 = CuboidalSurface( [0,0,0],[size,size,size])
# not supported yet
#s.addSurface( box1 )

#P = Species( 'P', 1e-12, 5e-8 )
P = Species( 'P', 1e-12, 3e-9 ) #hemo
s.addSpecies( P )

s.setAllRepulsive()

s.throwInParticles( 'P', 60*125, box1 )

l = Logger( s, 'simple' )
l.setParticleOutput( ('P',) )
l.setInterval( 1e-2 )
l.log()

#while s.t < 10:
#    s.step()
#    l.log()
    


def profrun():
    for i in range( 100 ):
        s.step()

import profile
profile.run('profrun()', 'fooprof')
import pstats
pstats.Stats('fooprof').sort_stats('time').print_stats(30)


sys.exit(1)


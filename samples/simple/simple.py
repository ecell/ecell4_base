#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

size = 1e-6
w = World(size, 3)
s = EGFRDSimulator(w)


box1 = CuboidalRegion([0,0,0],[size,size,size])
# not supported yet
#s.addSurface(box1)

#P = Species('P', 1e-12, 5e-8)
P = Species('P', 1e-12, 3e-9) #hemo
s.addSpecies(P)

s.setAllRepulsive()

s.throwInParticles(P, 60, box1)

l = Logger(s, 'simple')
l.setParticleOutput(('P', ))
l.setParticleOutInterval(3.33e-4)
l.log()

while s.t < .1:
    s.step()
    l.log()
    


# def profrun():
#     for i in range(100):
#         s.step()

# import profile
# profile.run('profrun()', 'fooprof')
# import pstats
# pstats.Stats('fooprof').sort_stats('time').print_stats(30)


# sys.exit(1)


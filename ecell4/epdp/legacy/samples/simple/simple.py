#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

size = 1e-6
w = World(size, 3)
s = EGFRDSimulator(w)


box1 = CuboidalRegion([0,0,0],[size,size,size])
# not supported yet
#s.add_surface(box1)

#P = Species('P', 1e-12, 5e-8)
P = Species('P', 1e-12, 3e-9) #hemo
s.add_species(P)

s.set_all_repulsive()

s.throw_in_particles(P, 60, box1)

l = Logger('simple')
interrupter = FixedIntervalInterrupter(s, 3.33e-4, l)

l.start(s)
while s.t < .1:
    interrupter.step()

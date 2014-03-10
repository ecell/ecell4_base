#!/usr/bin/env python

'''
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.3.out 1.25e-2 20000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.2.out 1.25e-3 20000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.1.out 1.25e-4 7000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.0.out 1.25e-5 5000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.-1.out 1.25e-6 2000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.-2.out 1.25e-7 2000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.-3.out 1.25e-8 1000000 &
'''

import sys
from egfrd import *
from bd import *
import model
import gfrdbase
import _gfrd
import myrandom

def run(outfilename, T, N):

    outfile = open(outfilename, 'w')

    for i in xrange(N):
        d, t = singlerun(T)
        outfile.write('%.18g\n' % d)
        outfile.flush()
        #print i
        #print d, t
        assert d == 0 or t == T


    outfile.close()



def singlerun(T):

    sigma = 5e-9
    r0 = sigma
    D = 1e-12
    D_tot = D * 2

    kf = 100 * sigma * D_tot

    m = model.ParticleModel(1e-3)

    A = model.Species('A', D, sigma/2)
    m.add_species_type(A)
    B = model.Species('B', D, sigma/2)
    m.add_species_type(B)
    C = model.Species('C', D, sigma/2)
    m.add_species_type(C)

    r1 = model.create_binding_reaction_rule(A, B, C, kf)
    m.network_rules.add_reaction_rule(r1)

    w = gfrdbase.create_world(m, 3)
    nrw = gfrdbase.create_network_rules_wrapper(m)
    s = _gfrd._EGFRDSimulator(w, nrw, myrandom.rng)

    class check_reactions:
        def __init__(self):
            self.reactions = []

        def __call__(self, ri):
            self.reactions.append(ri)

    cr = check_reactions()
    s.reaction_recorder = cr

    pid1 = gfrdbase.place_particle(w, A, [0,0,0])[0]
    pid2 = gfrdbase.place_particle(w, B, [float(A['radius']) + 
                                          float(B['radius'])+1e-23,0,0])[0]

    end_time = T

    while s.step(end_time):
        if len(cr.reactions) != 0:
            cr.reactions = []
            return 0, s.t

    p1 = s.world.get_particle(pid1)[1]
    p2 = s.world.get_particle(pid2)[1]

    distance = w.distance(p1.position, p2.position)

    return distance, s.t

    
if __name__ == '__main__':
    run(sys.argv[1], float(sys.argv[2]), int(sys.argv[3]))

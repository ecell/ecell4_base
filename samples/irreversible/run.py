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
import myrandom

def run(outfilename, T, N):

    outfile = open(outfilename, 'w')

    for i in range(N):
        d, t = singlerun2(T)
        outfile.write('%.18g\n' % d)
        outfile.flush()
        #print i
        #print d, t
        assert d == 0 or t == T


    outfile.close()



def singlerun1(T):
    #s.set_max_shell_size(1e-6)


    sigma = 5e-9
    r0 = sigma
    D = 2e-12
    D_tot = D
    kf = 10 * sigma * D_tot

    m = model.ParticleModel(1e-3)

    A = model.Species('A', 0.0, sigma/2)
    m.add_species_type(A)
    B = model.Species('B', D, sigma/2)
    m.add_species_type(B)
    C = model.Species('C', 0.0, sigma/2)
    m.add_species_type(C)

    r1 = model.create_binding_reaction_rule(A, B, C, kf)
    m.network_rules.add_reaction_rule(r1)

    w = gfrdbase.create_world(m, 3)
    nrw = gfrdbase.create_network_rules_wrapper(m)
    s = BDSimulator(w, myrandom.nrg, nrw)

    particleA = gfrdbase.place_particle(w, A, [0,0,0])
    particleB = gfrdbase.place_particle(w, B, [(float(A['radius']) + float(B['radius']))+1e-23,0,0])

    end_time = T
    s.step()

    while 1:
        next_time = s.get_next_time()
        if next_time > end_time:
            s.stop(end_time)
            break
        s.step()
        if s.last_reaction:
            print 'reaction'
            return 0.0, s.t

    distance = w.distance(s.get_position(particleB), s.get_position(particleA))

    return distance, s.t


def singlerun2(T):

    #s.set_user_max_shell_size(1e-7)
    #s.set_user_max_shell_size(1e-3)

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
    s = EGFRDSimulator(w, myrandom.rng, nrw)

    particleA = gfrdbase.place_particle(w, A, [0,0,0])
    particleB = gfrdbase.place_particle(w, B, [float(A['radius']) + float(B['radius'])+1e-23,0,0])

    end_time = T

    while 1:
        s.step()
        if s.last_reaction:
            #print 'reaction'
            return 0.0, s.t

        next_time = s.get_next_time()
        if next_time > end_time:
            s.stop(end_time)
            break

    distance = w.distance(s.get_position(particleA), s.get_position(particleB))

    return distance, s.t

    

if __name__ == '__main__':
    run(sys.argv[1], float(sys.argv[2]), int(sys.argv[3]))

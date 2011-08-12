#!/usr/bin/env python

'''
 PYTHONPATH=../.. python -O run.py single.0.out 5e-4 3e-8 100000
'''



from egfrd import *
import sys
import _gfrd
import model
import gfrdbase

import myrandom


def run(outfilename, T, S, N):
    print outfilename

    outfile = open(outfilename, 'w')

    for i in range(N):
        d, t = singlerun(T, S)
        outfile.write('%g\n' % d)

        #print i
        assert t == T

    outfile.close()



def singlerun(T, S):

    m = model.ParticleModel(1e-3)

    A = model.Species('A', 1e-12, 5e-9)
    m.add_species_type(A)

    w = gfrdbase.create_world(m, 3)
    nrw = gfrdbase.create_network_rules_wrapper(m)
    # s = EGFRDSimulator(w, myrandom.rng, nrw)
    # s.set_user_max_shell_size(S)
    s = _gfrd._EGFRDSimulator(w, nrw, myrandom.rng, 1, 1e-5, S)
    
    particleA = gfrdbase.place_particle(w, A, [0, 0, 0])

    end_time = T
    s.step()

    while 1:
        next_time = s.t + s.dt
        if next_time > end_time:
            # s.stop(end_time)
            s.step(end_time)
            break
        s.step()

    pos = w.get_particle(iter(w.get_particle_ids(A.id)).next())[1].position
    distance = w.distance([0,0,0], pos)
    return distance, s.t
    
def first(x):
    x = iter(x)
    try:
        return x.next()
    except StopIteration, e:
        return None

if __name__ == '__main__':
    run(sys.argv[1], float(sys.argv[2]), float(sys.argv[3]),
        int(sys.argv[4]))


#!/usr/bin/env python

'''
 PYTHONPATH=../.. python -O run.py single.0.out 5e-4 3e-8 100000
'''



from egfrd import *
import sys

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

    w = World(1e-3, 3)
    s = EGFRDSimulator(w)

    s.set_user_max_shell_size(S)

    m = ParticleModel()

    A = m.new_species_type('A', 1e-12, 5e-9)

    s.set_model(m)

    particleA = s.place_particle(A, [0,0,0])

    end_time = T
    s.step()

    while 1:
        next_time = s.get_next_time()
        if next_time > end_time:
            s.stop(end_time)
            break
        s.step()

    pos = s.get_position(A.id)
    distance = s.distance([0,0,0], pos)
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


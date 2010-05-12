#!/usr/bin/env python

from egfrd import *
from bd import *

from logger import *
import sys
import time


def run_single(T, V, N):

    print 'T =', T, '; V= ', V, '; N=', N
    

    L = math.pow(V * 1e-3, 1.0 / 3.0)

    matrix_size = max(3, int((3 * N) ** (1.0/3.0)))

    print 'matrix_size=', matrix_size
    
    w = World(L, matrix_size)
    s = EGFRDSimulator(w)
    #s = BDSimulator(w)

    box1 = CuboidalRegion([0,0,0],[L,L,L])

    D = 1e-12

    m = ParticleModel()

    A = m.new_species_type('A', D, 2.5e-9)
    m.set_all_repulsive()

    s.set_model(m)
    
    s.throw_in_particles(A, N, box1)
    print 'stir'

    t = 0
    stir_time = T * .1
    while 1:
        s.step()
        next_time = s.get_next_time()
        if next_time > stir_time:
            s.stop(stir_time)
            break

    print 'reset'
    s.reset()

    print 'run'
    run_time = T

    start = time.time()
    while s.t < run_time:
        s.step()
    end = time.time()

    timing = end - start

    steps = s.step_counter
    stepspersec = float(steps) / timing
    print 'steps (total)= ', steps
    print 'steps/sec= ', stepspersec, ', steps/N= ', float(steps) / N
    print 'TIMING:\n', timing, '\n'

    return end - start, steps, stepspersec


def run_single_bd(T, V, N, dt_factor):

    print 'T =', T, '; V= ', V, '; N=', N
    

    L = math.pow(V * 1e-3, 1.0 / 3.0)

    matrix_size = max(3, int((3 * N) ** (1.0/3.0)))

    print 'matrix_size=', matrix_size
    
    w = World(L, matrix_size)
    s = BDSimulator(w)

    box1 = CuboidalRegion([0,0,0],[L,L,L])

    D = 1e-12

    m = ParticleModel()

    A = m.new_species_type('A', D, 2.5e-9)
    m.set_all_repulsive()

    s.set_model(m)
    s.dt_factor = dt_factor
    
    s.throw_in_particles(A, N, box1)
    print 'stir'

    t = 0
    stir_time = T * .1
    while 1:
        s.step()
        next_time = s.get_next_time()
        if next_time > stir_time:
            s.stop(stir_time)
            break

    print 'reset'
    s.reset()

    print 'run'
    run_time = T

    start = time.time()
    while s.t < run_time:
        s.step()
    end = time.time()

    timing = end - start

    steps = s.step_counter
    stepspersec = float(steps) / timing
    print 'steps (total)= ', steps
    print 'steps/sec= ', stepspersec, ', steps/N= ', float(steps) / N
    print 'TIMING:\n', timing, '\n'

    return end - start, steps, stepspersec



if __name__ == '__main__':
    
    T = float(sys.argv[1])
    V = float(sys.argv[2])
    N = int(sys.argv[3])

    run_single(T, V, N)


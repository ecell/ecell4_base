#!/usr/bin/env python

from egfrd import *
from bd import *
from gfrdbase import *
from logger import *
import _gfrd
import model
import sys


N = 300

L = 5e-6
#L = 2e-6
#L = 5e-8
#L = 3e-7

m = model.ParticleModel(L)
S = model.Species('S', 1.5e-12, 5e-9)
P = model.Species('P', 1e-12, 7e-9)
m.add_species_type(S)
m.add_species_type(P)
r1 = model.create_binding_reaction_rule(S, S, P, 1e7 / N_A)
r2 = model.create_unbinding_reaction_rule(P, S, S, 1e3)
m.network_rules.add_reaction_rule(r1)
m.network_rules.add_reaction_rule(r2)
m.set_all_repulsive()

world = create_world(m, int((N * 6) ** (1. / 3.)))
nrw = _gfrd.NetworkRulesWrapper(m.network_rules)
s = EGFRDSimulator(world, myrandom.rng, nrw)
#s = BDSimulator(world. myrandom.rng, nrw)

throw_in_particles(s.world, S, N / 2)
throw_in_particles(s.world, P, N / 2)


#l = Logger('dimer')
l = None
interrupter = None

if l is not None:
    interrupter = FixedIntervalInterrupter(s, 1e-7, l.log)

import myrandom
myrandom.seed(0)


def profrun():
    if l is not None:
        l.start(s)
    for _ in xrange(15000):
        if interrupter is not None:
            interrupter.step()
        else:
            s.step()

PROFMODE = True

if PROFMODE:
    try:
        import cProfile as profile
    except:
        import profile
    profile.run('profrun()', 'fooprof')
    s.print_report()

    import pstats
    pstats.Stats('fooprof').sort_stats('time').print_stats(40)

else:
    profrun()
    s.print_report()

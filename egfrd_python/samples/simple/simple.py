import csv
import sys
import time

from ecell4.core import *
from ecell4.egfrd import *


def run():
    L = 1e-6
    volume = L * L * L
    N = 60

    D, radius = "1e-12", "2.5e-9"

    k2, U = 0.1, 0.5
    k1 = k2 * volume * (1 - U) / (U * U * N)

    sp1, sp2, sp3 = Species("A"), Species("B"), Species("C")
    species_list = [sp1, sp2, sp3]
    for sp in species_list:
        sp.set_attribute("D", D)
        sp.set_attribute("radius", radius)

    m = NetworkModel()
    m.add_species_attribute(sp1)
    m.add_species_attribute(sp2)
    m.add_species_attribute(sp3)
    m.add_reaction_rule(
        create_binding_reaction_rule(sp1, sp2, sp3, k1))
    m.add_reaction_rule(
        create_unbinding_reaction_rule(sp3, sp1, sp2, k2))

    rng = GSLRandomNumberGenerator()
    rng.seed(time.time())

    matrix_size = 3
    w = EGFRDWorld(L, matrix_size, rng)
    w.add_molecules(sp3, N)

    sim = EGFRDSimulatorWrapper(m, w)
    # sim.initialize()

    next_time, dt = 0.0, 0.02

    writer = csv.writer(sys.stdout, delimiter='\t')
    writer.writerow(['#t'] + [sp.name() for sp in species_list])
    writer.writerow(
        ['%.6e' % sim.t()]
        + ['%d' % w.num_molecules(sp) for sp in species_list])
    for i in xrange(10000):
        next_time += dt
        while sim.step(next_time):
            pass

        writer.writerow(
            ['%.6e' % sim.t()]
            + ['%d' % w.num_molecules(sp) for sp in species_list])


if __name__ == "__main__":
    run()

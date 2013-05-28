import csv
import sys

from ecell4.core import *
from ecell4.spatiocyte import *


def run():
    voxel_radius = 1e-8
    D = "1e-12"
    N = 60
    L = 1e-7
    volume = L * L * L
    k2, U = 0.1, 0.5
    k1 = k2 * volume * (1 - U) * (U * U * N)

    sp1, sp2, sp3 = Species("A"), Species("B"), Species("C")
    species_list = [sp1, sp2, sp3]
    for sp in species_list:
        sp.set_attribute("D", D)

    m = NetworkModel()
    m.add_species(sp1)
    m.add_species(sp2)
    m.add_species(sp3)
    m.add_reaction_rule(
        create_binding_reaction_rule(sp1, sp2, sp3, k1))
    m.add_reaction_rule(
        create_unbinding_reaction_rule(sp3, sp1, sp2, k2))

    w = SpatiocyteWorld(Position3(L, L, L), voxel_radius)
    w.add_molecules(sp3, N)

    sim = SpatiocyteSimulator(m, w)
    sim.initialize()

    next_time, dt = 0.0, 0.02

    writer = csv.writer(sys.stdout, delimiter='\t')
    writer.writerow(['#t'] + [sp.name() for sp in species_list])
    writer.writerow(
        ['%.6e' % sim.t()]
        + ['%d' % w.num_molecules(sp) for sp in species_list])
    for i in xrange(1000):
        next_time += dt
        while sim.step(next_time):
            pass

        writer.writerow(
            ['%.6e' % sim.t()]
            + ['%d' % w.num_molecules(sp) for sp in species_list])


if __name__ == "__main__":
    run()

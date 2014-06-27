import time

from ecell4.core import *
from ecell4.gillespie import *


def run():
    sp1, sp2, sp3 = Species("A"), Species("B"), Species("C")
    kf, kr = 0.25, 1.0

    rr1 = create_binding_reaction_rule(sp1, sp2, sp3, kf)
    rr2 = create_unbinding_reaction_rule(sp3, sp1, sp2, kr)

    m = NetworkModel()
    m.add_species_attribute(sp1)
    m.add_species_attribute(sp2)
    m.add_species_attribute(sp3)
    m.add_reaction_rule(rr1)
    m.add_reaction_rule(rr2)

    rng = GSLRandomNumberGenerator()
    rng.seed(int(time.time()))

    L = 1.0
    w = GillespieWorld(Position3(L, L, L), rng)
    # w.add_species(sp1)
    # w.add_species(sp2)
    # w.add_species(sp3)
    w.add_molecules(sp3, 10)
    w.save("test_gillespie.h5")

    sim = GillespieSimulator(m, w)

    print "t = %g, A = %g, B = %g, C = %g" % (
        sim.t(), w.num_molecules(sp1), w.num_molecules(sp2),
        w.num_molecules(sp3))
    for i in range(100):
        sim.step()
        print "t = %g, A = %g, B = %g, C = %g" % (
            sim.t(), w.num_molecules(sp1), w.num_molecules(sp2),
            w.num_molecules(sp3))


if __name__ == "__main__":
    run()

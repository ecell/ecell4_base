import time

from ecell4.core import *
from ecell4.gillespie import *


def run():
    sp1, sp2, sp3 = Species("A"), Species("B"), Species("C")
    kf, kr = 0.25, 1.0

    rr1 = ReactionRule()
    rr1.set_k(kf)
    rr1.add_reactant(sp1)
    rr1.add_reactant(sp2)
    rr1.add_product(sp3)

    rr2 = ReactionRule()
    rr2.set_k(kr)
    rr2.add_reactant(sp3)
    rr2.add_product(sp1)
    rr2.add_product(sp2)

    m = NetworkModel()
    m.add_species(sp1)
    m.add_species(sp2)
    m.add_species(sp3)
    m.add_reaction_rule(rr1)
    m.add_reaction_rule(rr2)

    rng = GSLRandomNumberGenerator()
    rng.seed(int(time.time()))

    volume = 1.0
    w = GillespieWorld(volume, rng)
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

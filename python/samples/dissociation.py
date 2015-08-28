from ecell4.core import *
from ecell4.ode import *


def singlerun():
    L, k1 = 1.0, 1.0
    sp1, sp2, sp3 = Species("A"), Species("B"), Species("C")
    rr1 = create_unbinding_reaction_rule(sp1, sp2, sp3, k1)

    m = NetworkModel()
    m.add_species_attribute(sp1)
    m.add_species_attribute(sp2)
    m.add_species_attribute(sp3)
    m.add_reaction_rule(rr1)

    w = ODEWorld(Real3(L, L, L))
    w.add_molecules(sp1, 60)

    target = ODESimulator(m, w)
    target.initialize()

    next_time = 0.0
    dt = 0.01

    print "t = %g\t A = %g\t B = %g\t C = %g" % (
        target.t(), w.num_molecules(sp1), w.num_molecules(sp2),
        w.num_molecules(sp3))
    for i in range(200):
        next_time += dt
        target.step(next_time)
        print "t = %g\t A = %g\t B = %g\t C = %g" % (
            target.t(), w.num_molecules(sp1), w.num_molecules(sp2),
            w.num_molecules(sp3))


if __name__ == "__main__":
    singlerun()

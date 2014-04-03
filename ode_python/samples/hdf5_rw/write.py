
from ecell4.core import *
from ecell4.ode import *

import h5py

def run():
    volume = 0.5
    sp1, sp2, sp3 = Species("A"), Species("B"), Species("C")
    k = 1.0
    rr1 = create_unbinding_reaction_rule(sp1, sp2, sp3, k)

    m = NetworkModel()
    m.add_species_attribute(sp1)
    m.add_species_attribute(sp2)
    m.add_species_attribute(sp3)
    m.add_reaction_rule(rr1)

    w = ODEWorld(volume)
    w.add_molecules(sp1, 60)
    w.save("test_dissociation.h5")

    target = ODESimulator(m, w)

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
        h5filename = "test_dissociation_{}.h5".format(target.t())
        if i == 99: # t == 1.0
            w.save(h5filename)

if __name__ == "__main__":
    run()


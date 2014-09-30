import numpy
import time
from ecell4.core import *

# from ecell4.lattice import (
#     LatticeWorld as world_type, LatticeSimulator as simulator_type)
from ecell4.gillespie import (
    GillespieWorld as world_type, GillespieSimulator as simulator_type)
# from ecell4.ode import (
#     ODEWorld as world_type, ODESimulator as simulator_type)


def singlerun(seed):
    L, voxel_radius = 1e-6, 2.5e-9

    radius, D = "2.5e-9", "1e-12"
    sp1 = Species("A", radius, D)
    sp2 = Species("B", radius, D)
    sp3 = Species("C", radius, D)

    N, kd, U = 60, 0.5, 0.5
    ka = kd * (L * L * L) * (1 - U) / (U * U * N)
    kon, koff = ka, kd
    # kD = 4 * numpy.pi * 4 * float(radius) * float(D)
    # kon = kD * ka / (kD + ka)
    # koff = kd * kon / ka
    rr1 = create_unbinding_reaction_rule(sp1, sp2, sp3, koff)
    rr2 = create_binding_reaction_rule(sp2, sp3, sp1, kon)

    m = NetworkModel()
    m.add_species_attribute(sp1)
    m.add_species_attribute(sp2)
    m.add_species_attribute(sp3)
    m.add_reaction_rule(rr1)
    m.add_reaction_rule(rr2)

    rng = GSLRandomNumberGenerator()
    rng.seed(seed)

    # w = world_type(Position3(L, L, L), voxel_radius, rng) # lattice
    w = world_type(Position3(L, L, L), rng) # gillespie
    # w = world_type(Position3(L, L, L)) # ode
    w.add_molecules(sp1, N)
    # w.save("test.h5")

    sim = simulator_type(m, w)
    sim.initialize()

    data = []

    def log(mode="a"):
        t, N1, N2, N3 = (sim.t(),
            w.num_molecules(sp1), w.num_molecules(sp2), w.num_molecules(sp3))
        print "%e\t%g\t%g\t%g" % (t, N1, N2, N3)
        data.append(numpy.array([t, N1, N2, N3]))

    next_time, dt = 0.0, 0.05
    log()
    for i in range(100):
        next_time += dt
        while sim.step(next_time):
            pass
        log()

    return numpy.asarray(data)

def run(filename="test.dat", num_trials=10):
    singlerun(0)

    # data = singlerun(0)
    # for i in range(1, num_trials):
    #     data += singlerun(i)
    # numpy.savetxt(filename, data / num_trials)


if __name__ == "__main__":
    import sys


    if len(sys.argv) > 2:
        run(sys.argv[1], int(sys.argv[2]))
    elif len(sys.argv) == 2:
        run(sys.argv[1])
    else:
        run()

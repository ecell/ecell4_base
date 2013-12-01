from ecell4.core import *
from ecell4.ode import *
from ecell4.load_ode_world import *

import h5py

def run():
    filename = "test_dissociation_1.0.h5"

    sp1, sp2, sp3 = Species("A"), Species("B"), Species("C")
    k = 1.0
    rr1 = create_unbinding_reaction_rule(sp1, sp2, sp3, k)
    m = NetworkModel()
    m.add_species_attribute(sp1)
    m.add_species_attribute(sp2)
    m.add_species_attribute(sp3)
    m.add_reaction_rule(rr1)

    w = hdf5load_ode(filename) 
    target = ODESimulator(m, w) 
    
    dt = 0.01
    next_time = w.t() + dt
    print w.volume()

    print "t = %g\t A = %g\t B = %g\t C = %g" % (
        target.t(), w.num_molecules(sp1), w.num_molecules(sp2),
        w.num_molecules(sp3))
    for i in range(200):
        next_time += dt
        target.step(next_time)
        print "t = %g\t A = %g\t B = %g\t C = %g" % (
            target.t(), w.num_molecules(sp1), w.num_molecules(sp2),
            w.num_molecules(sp3))

run()

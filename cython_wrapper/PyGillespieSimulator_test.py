import PyEcell4
import PyGillespie

def run():
    rng = PyEcell4.PyRandomNumberGenerator()
    sp1 = PyEcell4.PySpecies("A")
    sp2 = PyEcell4.PySpecies("B")

    rr1 = PyEcell4.PyReactionRule()
    rr1.set_k(5.001)
    rr1.add_reactant(sp1)
    rr1.add_product(sp2)

    m = PyEcell4.PyNetworkModel()
    m.add_reaction_rule(rr1)

    w = PyGillespie.PyGillespieWorld(1.0)
    w.add_species(sp1)
    w.add_species(sp2)

    w.add_molecules(sp1, 10)
    w.add_molecules(sp2, 10)
    m.add_species(sp1)
    m.add_species(sp2)

    sim = PyGillespie.PyGillespieSimulator(m, w, rng)
    #sim.save_hdf5_init("cython_wrapper_test.hdf5")

    for i in range(10):
        sim.step()
        print "t: ",  w.t() , "\tA: " , w.num_molecules(sp1) , "\tB: " , w.num_molecules(sp2) 
        #sim.save_hdf5()


run()

import ecell4
import Gillespie
import time

def run():
    rng = ecell4.RandomNumberGenerator()
    rng.seed(int(time.time()))
    sp1 = ecell4.Species("A")
    sp2 = ecell4.Species("B")

    rr1 = ecell4.ReactionRule()
    rr1.set_k(5.001)
    rr1.add_reactant(sp1)
    rr1.add_product(sp2)

    m = ecell4.NetworkModel()
    m.add_reaction_rule(rr1)

    w = Gillespie.GillespieWorld(1.0)
    w.add_species(sp1)
    w.add_species(sp2)

    w.add_molecules(sp1, 10)
    w.add_molecules(sp2, 10)
    m.add_species(sp1)
    m.add_species(sp2)

    sim = Gillespie.GillespieSimulator(m, w, rng)
    #sim.save_hdf5_init("cython_wrapper_test.hdf5")

    for i in range(10):
        sim.step()
        print "t: ",  w.t() , "\tA: " , w.num_molecules(sp1) , "\tB: " , w.num_molecules(sp2) 
        #sim.save_hdf5()


run()

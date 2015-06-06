from ecell4.core import *
from ecell4.ode import *

def singlerun():
    L = 1e-16
    edge_length = Real3(L, L, L)
    volume = L * L * L
    N = 60
    ka, U = 0.1, 0.5
    sp1, sp2, sp3 = Species("A"), Species("B"), Species("C")
    rr1 = ODEReactionRule()
    rr1.add_reactant(sp1, 1.0)
    rr1.add_product(sp2, 1.0)
    rr1.add_product(sp3, 1.0)
    rl1 = ODERatelawMassAction(ka)
    rr1.set_ratelaw( rl1 )
    #rr1.set_k(ka)
    print rr1.has_ratelaw()

    rr2 = ODEReactionRule()
    rr2.add_reactant(sp2, 1.0)
    rr2.add_reactant(sp3, 1.0)
    rr2.add_product(sp1, 1.0)
    kd = ka * volume * (1 - U) / (U * U * N)
    rl2 = ODERatelawMassAction(kd)
    rr2.set_ratelaw( rl2 )
    #rr2.set_k(kd)
    print rr2.has_ratelaw()

    m = ODENetworkModel()
    m.add_reaction_rule(rr1)
    m.add_reaction_rule(rr2)
    for r in m.ode_reaction_rules():
        print r.as_string()
    
    w = ODEWorld(edge_length)
    w.add_molecules(sp1, N)
    
    sim = ODESimulator2(m, w)
    next_time, dt = 0.0, 0.01
    
    print "{}\t{}\t{}\t{}".format(
        sim.t(), w.num_molecules(sp1), w.num_molecules(sp2), w.num_molecules(sp3))
    for i in xrange(1000):
        next_time += dt
        sim.step(next_time)
        print "{}\t{}\t{}\t{}".format(
            sim.t(), w.num_molecules(sp1), w.num_molecules(sp2), w.num_molecules(sp3))

singlerun()
print "done"

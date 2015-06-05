from ecell4.core import *
from ecell4.ode import *

def singlerun():
    L = 1e-16
    edge_length = Real3(L, L, L)
    volume = L * L * L
    N = 60
    ka = 0.1
    U = 0.5
    sp1, sp2, sp3 = Species("A"), Species("B"), Species("C")
    rr1 = ODEReactionRule()
    rr1.add_reactant(sp1, 1.0)
    rr1.add_product(sp2, 1.0)
    rr1.add_product(sp3, 1.0)
    rr1.set_ratelaw( ODERatelawMassAction(ka) )

    rr2 = ODEReactionRule()
    rr2.add_reactant(sp2, 1.0)
    rr2.add_reactant(sp3, 1.0)
    rr2.add_product(sp1, 1.0)
    kd = ka * volume * (1 - U) / (U * U * N)
    rr2.set_ratelaw( ODERatelawMassAction(kd) )

    m = ODENetworkModel()
    m.add_reaction_rule(rr1)
    m.add_reaction_rule(rr2)
    
    for r in m.ode_reaction_rules():
        print r.as_string()

    

singlerun()
print "done"

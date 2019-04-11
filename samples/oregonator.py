from ecell4.core import *
from ecell4.ode import *

import sys

def generate_ode_reaction(leftside, rightside, k = 0.0):
    reaction = ODEReactionRule()
    for (coeff, sp) in leftside:
        if isinstance(coeff, float) and isinstance(sp, Species):
            reaction.add_reactant(sp, coeff)
    for (coeff, sp) in rightside:
        if isinstance(coeff, float) and isinstance(sp, Species):
            reaction.add_product(sp, coeff)
    if (k != 0.0):
        reaction.set_k(k)
    return reaction

# Variable(adjust is needed)
init_A = 6.0e-2
init_B = 0.0
init_P = 0.0
init_X = 5.025e-11
init_Y = 3.0e-2
init_Z = 4.8e-8

f = 1.0

k_1 = 1.34
k_2 = 1.6e+9
k_3 = 8.0e+3
k_4 = 4.0e+7
k_5 = 1.0

show_interval_steps = 100

def singlerun():
    L = 1.0
    edge_length = Real3(L, L, L)
    volume = L * L * L

    # Species
    A, B = Species("A"), Species("B")     
    X, Y, Z = Species("X"), Species("Y"), Species("Z")
    P = Species("P")
    
    # A + Y -> X + P
    rr1 = generate_ode_reaction( [(1.0, A), (1.0, Y)], [(1.0, X), (1.0, P)] )
    rr1.set_k(k_1)
    # X + Y -> P
    rr2 = generate_ode_reaction( [(1.0, X), (1.0, Y)], [(2.0, P)] )
    rr2.set_k(k_2)
    # A + X -> 2X + 2Z
    rr3 = generate_ode_reaction( [(1.0, A), (1.0, X)], [(2.0, X), (1.0, Z)] )
    rr3.set_k(k_3)
    # X + X -> A + P
    rr4 = generate_ode_reaction( [(2.0, X)], [(1.0, A), (1.0, P)] )
    rr4.set_k(k_4)
    # Z -> fY
    rr5 = generate_ode_reaction( [(1.0, Z)], [(f, Y)] )
    rr5.set_k(k_5)
    
    m = ODENetworkModel()
    for r in [rr1, rr2, rr3, rr4, rr5]:
        m.add_reaction_rule(r)
    for r in m.ode_reaction_rules():
        print r.as_string()
    
    w = ODEWorld(edge_length)
    w.set_value(A, init_A)
    w.set_value(B, init_B)
    w.set_value(P, init_P)
    w.set_value(X, init_X)
    w.set_value(Y, init_Y)
    w.set_value(Z, init_Z)
    
    sim = ODESimulator(m, w)
    next_time, dt = 0.0, 0.01
    
    print "{:5f},{:e},{:e},{:e}".format(
        sim.t(), w.get_value(X), w.get_value(Y), w.get_value(Z))
    for i in xrange(200000):
        next_time += dt
        sim.step(next_time)
        if i % show_interval_steps == 0:
            #sys.stderr.write("{}\n".format(sim.t()) )    
            print "{:5f},{:e},{:e},{:e}".format(
                sim.t(), w.get_value(X), w.get_value(Y), w.get_value(Z))

singlerun()
print "done"

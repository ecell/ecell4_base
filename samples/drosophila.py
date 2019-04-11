from ecell4.core import *
from ecell4.ode import *
from ecell4.util.decorator import *


def ratelaw_func(reactants, products, volume, t, rr):
    #import ipdb; ipdb.set_trace()
    flux = float(0.1 * volume)
    for c in reactants:
        flux *= c/volume
    return flux

def ratelaw_func2(reactants, products, volume, t, rr):
    ka, U = 0.1, 0.5
    N = 60
    kd = ka * volume * (1 - U) / (U * U * N)
    flux = float(kd * volume)
    for c in reactants:
        flux *= c/volume
    return flux

def singlerun1():
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
    # rl1 = ODERatelawMassAction(ka)
    # rr1.set_ratelaw_massaction(rl1)
    # rl_cb = ODERatelawCallback(ratelaw_func)
    # rr1.set_ratelaw(rl_cb)
    # rl_cb = ODERatelawCallback(lambda r, p, V, t, rr: ka * r[0])
    # rr1.set_ratelaw(rl_cb)
    # rr1.set_ratelaw(ODERatelawCallback(lambda r, p, V, t, rr: ka * r[0]))

    rr2 = ODEReactionRule()
    rr2.add_reactant(sp2, 1.0)
    rr2.add_reactant(sp3, 1.0)
    rr2.add_product(sp1, 1.0)
    kd = ka * volume * (1 - U) / (U * U * N)
    # rl2 = ODERatelawMassAction(kd)
    # rr2.set_ratelaw(rl2)
    # rl_cb2 = ODERatelawCallback(ratelaw_func2)
    # rr2.set_ratelaw(rl_cb2)
    # rl_cb2 = ODERatelawCallback(lambda r, p, V, t, rr: kd / V * r[0] * r[1])
    # rr2.set_ratelaw(rl_cb2)
    # rr2.set_ratelaw(ODERatelawCallback(lambda r, p, V, t, rr: kd / V * r[0] * r[1]))

    m = ODENetworkModel()
    m.add_reaction_rule(rr1)
    m.add_reaction_rule(rr2)
    for r in m.ode_reaction_rules():
        print("# {}".format(r.as_string()))

    w = ODEWorld(edge_length)
    w.add_molecules(sp1, N)

    sim = ODESimulator(m, w)

    # obs = FixedIntervalNumberObserver(0.01, ["A", "B", "C"])
    # sim.run(10.0, obs)
    # print(obs.data())

    next_time, dt = 0.0, 0.01
    print("{}\t{}\t{}\t{}".format(
        sim.t(), w.get_value(sp1), w.get_value(sp2), w.get_value(sp3)))
    for i in xrange(1000):
        next_time += dt
        sim.step(next_time)
        print("{}\t{}\t{}\t{}".format(
            sim.t(), w.get_value(sp1), w.get_value(sp2), w.get_value(sp3)))

def singlerun2():
    L = 1e-16
    edge_length = Real3(L, L, L)
    volume = L * L * L
    N = 60
    ka, U = 0.1, 0.5
    kd = ka * volume * (1 - U) / (U * U * N)

    with reaction_rules():
        A == B + C | (lambda r, p, V, t, rr: ka * r[0],
                      lambda r, p, V, t, rr: kd / V * r[0] * r[1])

    m = get_model()

    w = ODEWorld(edge_length)
    w.add_molecules(Species("A"), N)

    sim = ODESimulator(m, w)
    obs = FixedIntervalNumberObserver(0.01, ["A", "B", "C"])
    sim.run(20.0, obs)

    for data in obs.data():
        print("{}\t{}\t{}\t{}".format(*data))

# singlerun1()
singlerun2()
print("# done")

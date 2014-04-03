from ecell4.core import *
from ecell4.ode import *

def run():
    volume = float(1e-18)
    N = float(60)
    ka, U = float(0.1), float(0.5)

    sp1, sp2, sp3 = Species("A"), Species("B"), Species("C")
    rr1 = ReactionRule()
    rr1.set_k(ka)
    rr1.add_reactant(sp1)
    rr1.add_product(sp2)
    rr1.add_product(sp3)

    kd = (float)(ka * volume * (1 - U) / (U * U * N))
    rr2 = ReactionRule()
    rr2.set_k(kd)
    rr2.add_reactant(sp2)
    rr2.add_reactant(sp3)
    rr2.add_product(sp1)

    model = NetworkModel()
    model.add_species_attribute(sp1)
    model.add_species_attribute(sp2)
    model.add_species_attribute(sp3)
    model.add_reaction_rule(rr1)
    model.add_reaction_rule(rr2)

    world = ODEWorld(volume)
    world.add_molecules(sp1, N)

    next_time = float(0.0)
    dt = float(0.01)
    target = ODESimulator(model, world)

    print "t = %g\t A = %g\t B = %g\t C = %g" % (
        target.t(), world.num_molecules(sp1),
        world.num_molecules(sp2),
        world.num_molecules(sp3))
    for i in range(1000):
        next_time += dt
        target.step(next_time)
        print "t = %g\t A = %g\t B = %g\t C = %g" % (
                target.t(),
                world.num_molecules(sp1),
                world.num_molecules(sp2),
                world.num_molecules(sp3))


if __name__ == "__main__":
    run()

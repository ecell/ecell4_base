from ecell4.reaction_reader.decorator import species_attributes, reaction_rules


@species_attributes
def attributes():
    K | {'N': '120'}
    KK | {'N': '30'}
    PP | {'N': '30'}

@reaction_rules
def reactions(kon1, koff1, kcat1, kon2, koff2, kcat2):
    (K + KK == K_KK | (kon1, koff1)
        > Kp + KK | kcat1
        == Kp_KK | (kon2, koff2)
        > Kpp + KK | kcat2)

    (Kpp + PP == Kpp_PP | (kon1, koff1)
        > Kp + PP | kcat1
        == Kp_PP | (kon2, koff2)
        > K + PP | kcat2)


if __name__ == "__main__":
    import sys
    import csv

    import ecell4.core
    import ecell4.gillespie


    m = ecell4.core.NetworkModel()

    for sp in attributes():
        m.add_species_attribute(sp)

    rules_list = reactions(
        4.483455086786913e-20, 1.35, 1.5,
        9.299017957780264e-20, 1.73, 15.0)
    for rr in rules_list:
        m.add_reaction_rule(rr)

    species_list = m.list_species()

    rng = ecell4.core.GSLRandomNumberGenerator()
    rng.seed(0)

    volume = 1e-18
    w = ecell4.gillespie.GillespieWorld(volume, rng)
    # for sp in attributes():
    for sp in species_list:
        sp = m.apply_species_attributes(sp)
        if sp.has_attribute("N"):
            w.add_molecules(sp, int(sp.get_attribute("N")))

    sim = ecell4.gillespie.GillespieSimulator(m, w)

    writer = csv.writer(sys.stdout, delimiter='\t')
    writer.writerow(['#t'] + [sp.name() for sp in species_list])
    writer.writerow(
        ['%.6e' % sim.t()]
        + ['%d' % w.num_molecules(sp) for sp in species_list])
    for i in xrange(1000):
        sim.step()
        writer.writerow(
            ['%.6e' % sim.t()]
            + ['%d' % w.num_molecules(sp) for sp in species_list])

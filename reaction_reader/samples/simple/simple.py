from ecell4.reaction_reader.decorator import species_attributes_with_keys, reaction_rules


@species_attributes_with_keys('N')
def species():
    K | '120'
    KK | '30'
    PP | '30'

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

    species_list = species()
    for sp in species_list:
        m.add_species_attribute(sp)

    rules = reactions(
        4.483455086786913e-20, 1.35, 1.5,
        9.299017957780264e-20, 1.73, 15.0)
    for rr in rules:
        m.add_reaction_rule(rr)

    rng = ecell4.core.GSLRandomNumberGenerator()
    rng.seed(0)

    volume = 1e-18
    w = ecell4.gillespie.GillespieWorld(volume, rng)
    for sp in species_list:
        if sp.has_attribute("N"):
            w.add_molecules(sp, int(sp.get_attribute("N")))

    sim = ecell4.gillespie.GillespieSimulator(m, w)

    writer = csv.writer(sys.stdout, delimiter='\t')
    writer.writerow(['#t'] + [sp.name() for sp in m.list_species()])
    writer.writerow(
        ['%.6e' % sim.t()]
        + ['%d' % w.num_molecules(sp) for sp in m.list_species()])
    for i in xrange(1000):
        sim.step()
        writer.writerow(
            ['%.6e' % sim.t()]
            + ['%d' % w.num_molecules(sp) for sp in m.list_species()])

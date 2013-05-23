from ecell4.reaction_reader.decorator import species_attributes, reaction_rules


@species_attributes
def species():
    K | {'N': '120'}
    KK | {'N': '30'}
    PP | {'N': '30'}
    Kp | {}
    Kpp | {}
    K_KK | {}
    Kp_KK | {}
    Kpp_PP | {}
    Kp_PP | {}

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
    import ecell4.core
    import ecell4.gillespie


    m = ecell4.core.NetworkModel()

    species_list = species()
    for sp in species_list:
        m.add_species(sp)

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
            N = int(sp.get_attribute("N"))
            if N > 0:
                w.add_molecules(sp, N)

    sim = ecell4.gillespie.GillespieSimulator(m, w)

    derim = '\t'
    print '#t%s%s' % (derim, derim.join((sp.name() for sp in species_list)))
    print '%.6e%s%s' % (
        sim.t(), derim,
        derim.join(("%d" % w.num_molecules(sp) for sp in species_list)))
    for i in xrange(1000):
        sim.step()
        print '%.6e%s%s' % (
            sim.t(), derim,
            derim.join(("%d" % w.num_molecules(sp) for sp in species_list)))

from ecell4.reaction_reader.decorator import species_attributes, reaction_rules


@species_attributes
def species():
    K | {'N': '10'}
    KK | {'N': '20'}
    PP | {'N': '10'}

@reaction_rules
def reactions(kon, koff, kcat):
    (K + KK == K_KK | (kon, koff)
        > Kp + KK | kcat
        == Kp_KK | (kon, koff)
        > Kpp + KK | kcat)

    (Kpp + PP == Kpp_PP | (kon, koff)
        > Kp + PP | kcat
        == Kp_PP | (kon, koff)
        > K + PP | kcat)

@reaction_rules
def synthesis(kf, kr):
    ~A > A | kf
    A > ~A | kr


if __name__ == "__main__":
    rules = reactions(1, 2, 3)
    # rules = synthesis(1, 0.1)

    for i, rr in enumerate(rules):
        # print i + 1, rr
        print i + 1, (
            [sp.name() for sp in rr.reactants()],
            [sp.name() for sp in rr.products()],
            rr.k())

    # import ecell4.core
    # import ecell4.gillespie


    # m = ecell4.core.NetworkModel()

    # for sp in species():
    #     print sp.name(), int(sp.get_attribute("N")) if sp.has_attribute("N") else 0

    # species_list = []
    # for rr in rules:
    #     m.add_reaction_rule(rr)
    #     species_list.extend(rr.reactants())
    #     species_list.extend(rr.products())

    # species_list = sorted(set(species_list))
    # for sp in species_list:
    #     m.add_species(sp)

    # rng = ecell4.core.GSLRandomNumberGenerator()
    # rng.seed(0)

    # volume = 1.0
    # w = ecell4.gillespie.GillespieWorld(volume, rng)
    # # for sp in species_list:
    # #     w.add_molecules(sp, 10)

    # sim = ecell4.gillespie.GillespieSimulator(m, w)

    # derim = '\t'
    # print '#t%s%s' % (derim, derim.join((sp.name() for sp in species_list)))
    # print '%.6e%s%s' % (
    #     sim.t(), derim,
    #     derim.join(("%d" % w.num_molecules(sp) for sp in species_list)))
    # for i in xrange(100):
    #     sim.step()
    #     print '%.6e%s%s' % (
    #         sim.t(), derim,
    #         derim.join(("%d" % w.num_molecules(sp) for sp in species_list)))

from ecell4.reaction_reader.decorator import reaction_rules


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


if __name__ == "__main__":
    rules = reactions(1, 2, 3)

    species_list = []
    for i, rr in enumerate(rules):
        print i + 1, rr

    #     reactants, products = rr.reactants(), rr.products()
    #     species_list.extend(reactants)
    #     species_list.extend(products)

    #     print i + 1, (
    #         [sp.name() for sp in reactants],
    #         [sp.name() for sp in products],
    #         rr.k())

    # species_list = sorted(set(species_list))
    # print tuple(sp.name() for sp in species_list)

from ecell4.reaction_reader.decorator import reaction_rules


@reaction_rules
def reactions(kon, koff, kcat):
    (A + B == C | (kon, koff)
        > D | kcat)

    D > C | kcat


if __name__ == "__main__":
    rules = reactions(1, 2, 3)

    # species = []
    # for rr in rules:
    #     species.extend(rr.reactants())
    #     species.extend(rr.products())

    # print [sp.name() for sp in species]
    # species = sorted(set(species))
    # print [sp.name() for sp in species]

    for i, rr in enumerate(rules):
        print i + 1, rr
        # reactants, products = rr.reactants(), rr.products()
        # print (i + 1,
        #     [sp.name() for sp in reactants],
        #     [sp.name() for sp in products],
        #     rr.k())

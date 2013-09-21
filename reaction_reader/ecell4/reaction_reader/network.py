import species
import ecell4.core as core


def generate_NetworkModel(seeds, rules):
    model = core.NetworkModel()

    for sp in seeds:
        model.add_species_attribute(core.Species(str(sp)))

    for r_tuple in rules:
        rr = core.ReactionRule()
        for reactant in r_tuple[0]:
            rr.add_reactant(core.Species(str(reactant)))
        for product in r_tuple[1]:
            rr.add_product(core.Species(str(product)))

        # kinetic parameter
        for opt in r_tuple[2]:
            if (not isinstance(opt, bool)
                and isinstance(opt, (int, long, float, complex))):
                rr.set_k(opt)
                break
        else:
            raise RuntimeError, "No kinetic rate is defined. [%s]" % str(r_tuple)

        model.add_reaction_rule(rr)

    return model


if __name__ == '__main__':
    pass

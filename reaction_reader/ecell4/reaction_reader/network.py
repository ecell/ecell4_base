import itertools
import copy

import species


def check_stoichiometry(sp, max_stoich):
    for pttrn, num_subunits in max_stoich.items():
        if sp.count_subunits(pttrn) > num_subunits:
            return False
    return True

def dump_reaction(reactants, products):
    # reactants, products = reaction
    for sp in itertools.chain(reactants, products):
        sp.sort()

    retval = "+".join(sorted([str(sp) for sp in reactants]))
    retval += ">"
    retval += "+".join(sorted([str(sp) for sp in products]))
    return retval

def reaction_rule_match_recurse(
    rr, idx, seeds1, seeds2, reactants, contexts, ignore):
    if idx >= rr.num_reactants():
        pttrns = rr.generate_with_contexts(reactants, contexts)
        return [
            (copy.deepcopy(reactants), products, rr.options())
            for products in pttrns]
    elif idx == ignore[0]:
        return reaction_rule_match_recurse(
            rr, idx + 1, seeds1, seeds2, reactants + [ignore[1]],
            contexts, ignore)

    if idx < ignore[0]:
        seeds = itertools.chain(seeds1, seeds2)
    else:
        seeds = seeds2

    retval = []
    for sp in seeds:
        newcontexts = rr.match_partial(idx, sp, copy.deepcopy(contexts))
        if newcontexts is None or len(contexts) == 0:
            continue
        tmp = reaction_rule_match_recurse(
            rr, idx + 1, seeds1, seeds2, reactants + [sp], newcontexts, ignore)
        retval.extend(tmp)
    return retval

def generate_recurse2(seeds1, rules, seeds2, max_stoich={}):
    seeds = list(itertools.chain(seeds1, seeds2))
    newseeds, newreactions = [], []
    for rr in rules:
        num_reactants = rr.num_reactants()
        if num_reactants == 0:
            continue # skip rules for synthesis

        for i in range(num_reactants):
            for sp in seeds1:
                contexts = rr.match_partial(i, sp)
                if contexts is None or len(contexts) == 0:
                    continue

                reactions = reaction_rule_match_recurse(
                    rr, 0, seeds1, seeds2, [], contexts, ignore=(i, sp))

                newreactions.extend(reactions)
                for (reactants, products, opts) in reactions:
                    for newsp in products:
                        if (newsp not in seeds and newsp not in newseeds
                            and check_stoichiometry(newsp, max_stoich)):
                            newsp.sort()
                            newseeds.append(newsp)
    return (newseeds, seeds, newreactions)

def generate_recurse(seeds1, rules, seeds2, max_stoich):
    seeds = list(itertools.chain(seeds1, seeds2))
    newseeds, newreactions = [], []
    for sp1 in seeds1:
        for rr in rules:
            if rr.num_reactants() == 1:
                pttrns = rr.generate(sp1)
                # try:
                #     pttrns = rr.generate(sp1)
                # except Exception, e:
                #     print rr, sp1
                #     raise e
                if pttrns is not None and len(pttrns) > 0:
                    for products in pttrns:
                        newreactions.append(((sp1, ), products, rr.options()))

                    for newsp in itertools.chain(*pttrns):
                        if (newsp not in seeds and newsp not in newseeds
                            and check_stoichiometry(newsp, max_stoich)):
                            newsp.sort()
                            newseeds.append(newsp)
        for sp2 in seeds:
            for rr in rules:
                if rr.num_reactants() == 2:
                    pttrns = rr.generate(sp1, sp2)
                    # try:
                    #     pttrns = rr.generate(sp1, sp2)
                    # except Exception, e:
                    #     print rr, sp1, sp2
                    #     raise e
                    if pttrns is not None and len(pttrns) > 0:
                        for products in pttrns:
                            newreactions.append(((sp1, sp2), products, rr.options()))

                        for newsp in itertools.chain(*pttrns):
                            if (newsp not in seeds and newsp not in newseeds
                                and check_stoichiometry(newsp, max_stoich)):
                                newsp.sort()
                                newseeds.append(newsp)
        for sp2 in seeds2:
            for rr in rules:
                if rr.num_reactants() == 2:
                    pttrns = rr.generate(sp2, sp1)
                    # try:
                    #     pttrns = rr.generate(sp2, sp1)
                    # except Exception, e:
                    #     print rr, sp1, sp2
                    #     raise e
                    if pttrns is not None and len(pttrns) > 0:
                        for products in pttrns:
                            newreactions.append(((sp1, sp2), products, rr.options()))

                        for newsp in itertools.chain(*pttrns):
                            if (newsp not in seeds and newsp not in newseeds
                                and check_stoichiometry(newsp, max_stoich)):
                                newsp.sort()
                                newseeds.append(newsp)
    return (newseeds, seeds, newreactions)

def generate_reactions(newseeds, rules, max_iter=10, max_stoich={}):
    seeds, cnt, reactions = [], 0, []

    for rr in rules:
        if rr.num_reactants() == 0:
            reactions.append((rr.reactants(), rr.products(), rr.options()))
            for newsp in rr.products():
                if (newsp not in newseeds
                    and check_stoichiometry(newsp, max_stoich)):
                    newsp.sort()
                    newseeds.append(newsp)

    while len(newseeds) != 0 and cnt < max_iter:
        # print "[RESULT%d] %d seeds, %d newseeds, %d reactions." % (
        #     cnt, len(seeds), len(newseeds), len(reactions))
        # newseeds, seeds, newreactions = generate_recurse(
        newseeds, seeds, newreactions = generate_recurse2(
            newseeds, rules, seeds, max_stoich)
        reactions.extend(newreactions)
        cnt += 1
    # print "[RESULT%d] %d seeds, %d newseeds, %d reactions." % (
    #     cnt, len(seeds), len(newseeds), len(reactions))
    # print ""

    seeds.sort(key=str)
    # for i, sp in enumerate(seeds):
    #     print "%5d %s" % (i + 1, str(sp))
    # print ""

    # reactions = list(set([dump_reaction(reaction) for reaction in reactions]))
    dump_rrobj_map = dict()
    for r in reactions:
        s = dump_reaction(r[0], r[1])
        dump_rrobj_map[s] = r
    reactions = dump_rrobj_map.values()
    # for i, reaction in enumerate(reactions):
    #     print "%5d %s" % (i + 1, reaction)

    return seeds + newseeds, reactions

def generate_NetworkModel(seeds, rules):
    import ecell4.core as core


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

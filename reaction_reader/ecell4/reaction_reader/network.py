import species
import ecell4.core as core

def generateNetworkModel(seeds, rules):
    model = core.NetworkModel()
    
    for sp in seeds:
        model.add_species_attribute( core.Species(str(sp)) )
    for r_tuple in rules:
        rr = core.ReactionRule()
        for reactant in r_tuple[0]:
            rr.add_reactant(core.Species( str(reactant) ))
        for product in r_tuple[1]:
            rr.add_product(core.Species( str(product) ))
        opts = r_tuple[2]
        # kinetic parameter
        if isinstance(opts, list) and 0 < len(opts):
            rr.set_k(opts[0])
        else:
            rr.set_k(opts)
        model.add_reaction_rule(rr)
    return model

        

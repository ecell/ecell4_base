
from ecell4.reaction_reader.decorator2 import species_attributes, reaction_rules
from ecell4.reaction_reader.species import generate_reactions
from ecell4.reaction_reader.sbml_exporter import convert2SBML
import ecell4.core as core

import ecell4.ode as ode
from ecell4.reaction_reader.network import generate_NetworkModel
from egfr import attributegen, rulegen
newseeds = []
attrs = attributegen()
for i, (sp, attr) in enumerate(attrs):
    #print i, sp, attr
    newseeds.append(sp)
    #print ''
reaction_rules = rulegen()

seeds, rules = generate_reactions(newseeds , reaction_rules, max_iter = 5)

m = generate_NetworkModel(seeds, rules)
convert2SBML(m, attrs, "egfr_n5.xml")

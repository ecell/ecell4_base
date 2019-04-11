import ecell4.core as core
import ecell4.ode as ode

from ecell4.reaction_reader.decorator2 import species_attributes, reaction_rules
from ecell4.reaction_reader.network import generate_reactions, generate_NetworkModel

from egfr import attributegen, rulegen


newseeds = []
attrs = attributegen()
for i, (sp, attr) in enumerate(attrs):
    #print i, sp, attr
    newseeds.append(sp)
    #print ''
reaction_rules = rulegen()

seeds, rules = generate_reactions(newseeds, reaction_rules, max_iter=3)
m = generate_NetworkModel(seeds, rules)
w = ode.ODEWorld(1.0)
for (sp, attr) in attrs:
    w.add_molecules(core.Species(str(sp)), attr)

target = ode.ODESimulator(m, w)
next_time = 0.0
dt = 0.01

for i in range(100):
    next_time += dt
    print "{}\t{} = {}".format(target.t(), str(seeds[0]), w.num_molecules(core.Species(str(seeds[0]))))
    target.step(next_time)

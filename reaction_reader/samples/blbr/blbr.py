from ecell4.reaction_reader.decorator2 import species_attributes, reaction_rules
from ecell4.reaction_reader.species import generate_reactions

@species_attributes
def attributegen():
    R(r,r) | R0
    L(l,l) | L0

@reaction_rules
def rulegen():
    # Ligand addition
    R(r) + L(l,l) == R(r^1).L(l^1,l) | (kp1,km1)

    # Chain elongation
    R(r) + L(l,l^_) == R(r^1).L(l^1,l^_) | (kp2,km2)

    # Ring closure
    R(r).L(l) == R(r^1).L(l^1) | (kp3, km3)

if __name__ == "__main__":
    newseeds = []
    for i, (sp, attr) in enumerate(attributegen()):
        print i, sp, attr
        newseeds.append(sp)
    print ''

    rules = rulegen()
    for i, rr in enumerate(rules):
        print i, rr
    print ''

    generate_reactions(newseeds, rules)

# setOption("SpeciesLabel","HNauty")
# begin model
# begin parameters
#     kp1  1
#     km1  1
#     kp2  1
#     km2  1
#     kp3  1
#     km3  1
#     R0   3e5
#     L0   3e5
# end parameters
# 
# begin seed species
#     R(r,r) R0
#     L(l,l) L0
# end seed species
# 
# begin reaction rules
#     # Ligand addition
#     R(r) + L(l,l) <-> R(r!1).L(l!1,l) kp1,km1
# 
#     # Chain elongation
#     R(r) + L(l,l!+) <-> R(r!1).L(l!1,l!+) kp2,km2
# 
#     # Ring closure
#     R(r).L(l) <-> R(r!1).L(l!1) kp3,km3
# end reaction rules
# end model
# 
# ## actions ##
# generate_network({overwrite=>1,max_stoich=>{R=>5,L=>5}})

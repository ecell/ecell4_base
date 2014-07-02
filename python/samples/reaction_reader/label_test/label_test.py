from ecell4.reaction_reader.decorator2 import species_attributes, reaction_rules
from ecell4.reaction_reader.species import generate_reactions 

from ecell4.reaction_reader.bng_exporter import check_label_containing_reaction , Convert2BNGManager

@species_attributes
def attributegen():
    R0 = 14
    L0 = 15
    #A(r1,r2) | R0
    #B(l1,l2) | L0
    A(x=P, phos=YT^1).C(bs^1) | R0 
    B(y=S, phos=pT^1).D(bs^1) | L0 


Kp1 = 1.0
Kp2 = 1.0

@reaction_rules
def rulegen():
    #A(x=_1) > B(y=_1) | Kp1
    A(x=_1) + B(y=_1) > A(x=_1^1).B(y=_1^1) | Kp2 
    A(x=_1, phos=_2) + B(y=_2) > A(x=_1, phos=_2^1).B(y=_1^1) | Kp2


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
    s = Convert2BNGManager(attributegen(), rulegen() )

    with open("expanded.bngl", "w") as fd:
        s.write_section_molecule_types(fd)
        s.write_section_seed_species(fd)
        s.write_section_reaction_rules(fd)
    #generate_reactions(newseeds, rules)

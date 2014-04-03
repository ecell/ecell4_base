from ecell4.reaction_reader.decorator2 import species_attributes, reaction_rules
from ecell4.reaction_reader.network import generate_reactions


@species_attributes
def attributegen():
    A(a) | 1
    B(b) | 1

@reaction_rules
def rulegen():
    A(a) + B(b) > A(a^1).B(b^1) | kp1
 
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


#begin model
#begin parameters
#    kp1  1
#    km1  100
#end parameters
#
#begin molecule types
#    A(a)
#    B(b)
#end molecule types 
#
#begin seed species
#    $A(a)  1
#    B(b)   1
#end seed species
#
#begin reaction rules
#    A(a) + B(b) -> A(a!1).B(b!1)  kp1
#end reaction rules
#end model
#
### actions ##
#generate_network({overwrite=>1})
#simulate({method=>"ode",t_end=>10,n_steps=>20})

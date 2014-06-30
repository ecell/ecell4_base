
from ecell4.reaction_reader.decorator2 import species_attributes, reaction_rules
from ecell4.reaction_reader.species import generate_reactions
from ecell4.reaction_reader.bng_exporter import Convert2BNGManager
import sys

@species_attributes
def attributegen():
    X(b=off) | 3
    Y(a=pY) | 4

@reaction_rules
def rulegen():
    X(a^_,b=off) > X(a^_,b=on) | 1 # X(a!+,b~off)
    #X(a^_1) + ...  
    X(a=_1,b=_2) > Y(a=_1) + Z(b=_2)  | 2# X(a%1,b%2)
    Y(a=_1) + Z(b=_1) > X(a=_1,b=_1) | 3
    _1(b=off) > _1(b=on) | 4
    _1(a) + _2(a) > _1(a^1)._2(a^1) | 5
    _1(a) + _1(a) > _1(a^1)._1(a^1) | 6
    #_1(a=_2) + ...

bng_mng = Convert2BNGManager(attributegen(), rulegen() )
bng_mng.write_section_molecule_types(sys.stdout)
bng_mng.write_section_seed_species(sys.stdout)
bng_mng.write_section_reaction_rules(sys.stdout)

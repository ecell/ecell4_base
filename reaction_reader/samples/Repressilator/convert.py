from ecell4.reaction_reader.decorator2 import species_attributes, reaction_rules
from ecell4.reaction_reader.species import generate_reactions 
from ecell4.reaction_reader.bng_exporter import Convert2BNGManager

from Repressilator import attributegen, rulegen
import math

rules = rulegen(
        6.022e23 ,1.4e-15, 1e9, 224, 9, 0.5, 5e-4, 0.167,
        math.log(2) / 120, math.log(2) / 600, 1e-4, 1000, 1000)

with open("export.bngl", "w") as f:
    bng_mng = Convert2BNGManager(attributegen(), rules)
    bng_mng.write_section_molecule_types(f)
    bng_mng.write_section_seed_species(f)
    bng_mng.write_section_reaction_rules(f)

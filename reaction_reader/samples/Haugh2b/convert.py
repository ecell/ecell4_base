from ecell4.reaction_reader.decorator2 import species_attributes, reaction_rules
from ecell4.reaction_reader.species import generate_reactions 
from ecell4.reaction_reader.bng_exporter import Convert2BNGManager

from Haugh2b import attributegen, rulegen


with open("new_export.bngl", "w") as f:
    rules = rulegen(1, 0.1, 1, 0.001, 0.1, 90, 10, 10, 99, 1, 100)
    bng_mng = Convert2BNGManager(attributegen(1), rules)
    bng_mng.write_section_molecule_types(f)
    bng_mng.write_section_seed_species(f)
    bng_mng.write_section_reaction_rules(f)

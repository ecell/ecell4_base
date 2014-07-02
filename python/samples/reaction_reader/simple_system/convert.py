from ecell4.reaction_reader.decorator2 import species_attributes, reaction_rules
from ecell4.reaction_reader.species import generate_reactions 
from ecell4.reaction_reader.bng_exporter import Convert2BNGManager

from simple_system import attributegen, rulegen

'''
with open("too_old_export.bngl", "w") as fd:
        export_bng(fd, attributegen(), rulegen())
'''

#with open("old_export.bngl", "w") as fd:
# convert2bng_moleculetypes(fd, rulegen() )
# convert2bng_seed_species(fd, attributegen() )
# convert2bng_reaction_rules(fd, rulegen() )

with open("new_export.bngl", "w") as f:
    bng_mng = Convert2BNGManager(attributegen(), rulegen())
    bng_mng.write_section_molecule_types(f)
    bng_mng.write_section_seed_species(f)
    bng_mng.write_section_reaction_rules(f)

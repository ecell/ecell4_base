from ecell4.reaction_reader.decorator2 import species_attributes, reaction_rules
from ecell4.reaction_reader.species import generate_reactions 
from ecell4.reaction_reader.bng_exporter import Convert2BNGManager

from SHP2_base_model import attributegen, rulegen

'''
with open("too_old_export.bngl", "w") as fd:
        export_bng(fd, attributegen(), rulegen())
'''

#with open("old_export.bngl", "w") as fd:
# convert2bng_moleculetypes(fd, rulegen() )
# convert2bng_seed_species(fd, attributegen() )
# convert2bng_reaction_rules(fd, rulegen() )

with open("new_export.bngl", "w") as f:
    bng_mng = Convert2BNGManager(attributegen(), rulegen(1000, 10, 500, 1, 1, 1, 1, 0.1, 1, 10, 1, 1000, 100, 1000, 1000, 100, 100, 100, 1000, 100, 100, 1000, 0.025, 0.05))
    bng_mng.write_section_molecule_types(f)
    bng_mng.write_section_seed_species(f)
    bng_mng.write_section_reaction_rules(f)

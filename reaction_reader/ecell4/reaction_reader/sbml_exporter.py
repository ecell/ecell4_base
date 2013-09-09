
import sys
import os.path
from libsbml import *
import ecell4.core as core
import species 

Level = 2
Version = 4

def convert2SBML(nw_model, sp_attrs):
    level = 2
    version = 4
    sbmlDoc = SBMLDocument(level, version)
    if not isinstance(nw_model, core.NetworkModel):
        return
    
    # Model {{{
    model = sbmlDoc.createModel()
    model.setId("XXX")  # XXX
    #   }}}

    # Define Units   {{{
    unitdef = model.createUnitDefinition()
    unitdef.setId("per_second")
    # }}}

    # Compartment {{{
    compName = "cytosol"
    comp = model.createCompartment()
    comp.setId(compName)

    comp.setSize(1e-14)    # XXX
    # }}}

    # Draw Species {{{
    # first, seed species
    all_species = dict()

    sp_index = 0
    for sp in nw_model.list_species():
        all_species[ sp.name() ] = (sp_index, 0)
        sp_index += 1
    for (sp, attr) in sp_attrs:
        if not all_species.has_key(sp):
            all_species[ str(sp) ] = (sp_index, attr)
            sp_index += 1

    for (sp, attr) in all_species.items():
        sp = model.createSpecies()
        sp.setId("sp{}".format(sp_index))
        sp.setName(str(sp))
    import ipdb; ipdb.set_trace()


    # Draw Reactions {{{
    #   reaction = model.createReaction()
    #   reaction.setId("r{}".format(r_id) )

    r_id = 0
    for rr in model.query_reaction_rules():
        reaction = model.createReaction()
        reaction.setId("r{}".format(r_index))
        r_id += 1

    # Kinetic Law 
    k1 = reaction.createKineticLaw()

    #   }}}
    return sbmlDoc


#sbmlDoc = convert2SBML()
#SBMLok = sbmlDoc.checkInternalConsistency();
'''
if SBMLok == 0:
    print "success"
else:
    print "fail"
'''

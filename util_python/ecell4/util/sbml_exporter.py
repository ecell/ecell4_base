
import sys
import os.path
from libsbml import *
import ecell4.core as core
import ecell4.reaction_reader.species 

Level = 2
Version = 4

def convert2SBML(nw_model, sp_attrs, fname):
    level = 2
    version = 4
    sbmlDoc = SBMLDocument(level, version)
    if not isinstance(nw_model, core.NetworkModel):
        return
    
    # Model {{{
    model = sbmlDoc.createModel()
    model.setId("E-Cell4_Model")  # XXX
    #   }}}

    # Define Units   {{{
    unitdef = model.createUnitDefinition()
    unitdef.setId("per_second")
    # }}}

    # Compartment {{{
    compName = "cytosol"
    comp = model.createCompartment()
    comp.setId(compName)
     
    comp.setSize(1.0)    # XXX
    # }}}

    # Draw Species {{{
    # first, seed species
    all_species = dict()

    sp_index = 0
    # Distribute IDs for each species
    for (sp, attr) in sp_attrs:
        all_species[ str(sp) ] = ("sp{}".format(sp_index), attr)
        sp_index += 1

    for sp in nw_model.list_species():
        if not all_species.has_key( sp.name() ):
            all_species[ sp.name() ] = ("sp{}".format(sp_index), 0)
            sp_index += 1
    
    # Add to SBML model object

    for (sp, (sp_id, attr) ) in all_species.items():
        sbml_sp = model.createSpecies()
        sbml_sp.setId(sp_id)
        sbml_sp.setName(str(sp))
        if isinstance(attr, (int, long, float, complex)):
            sbml_sp.setInitialAmount(attr)
        sbml_sp.setCompartment(compName)

    # Draw Reactions {{{
    r_index = 0
    for rr in nw_model.reaction_rules():
        sbml_reaction = model.createReaction()
        sbml_reaction.setId("r{}".format(r_index))
        for reactant in rr.reactants():
            sbml_spr = sbml_reaction.createReactant()
            (sbml_spid, attr) = all_species[ reactant.name()]
            #sbml_spr.setSpecies( all_species[ reactant.name() ][0] )
            sbml_spr.setSpecies( sbml_spid)
        for product in rr.products():
            sbml_spr = sbml_reaction.createProduct()
            (sbml_spid, attr) = all_species[ product.name() ]
            sbml_spr.setSpecies( sbml_spid)
        r_index += 1
        # Kinetic Law 
        #k1 = sbml_reaction.createKineticLaw()

    #   }}}
    writeSBML(sbmlDoc, fname)
    return


#sbmlDoc = convert2SBML()
#SBMLok = sbmlDoc.checkInternalConsistency();
'''
if SBMLok == 0:
    print "success"
else:
    print "fail"
'''

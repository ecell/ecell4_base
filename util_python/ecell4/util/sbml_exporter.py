
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
    #import ipdb; ipdb.set_trace()
    for (sp, attr) in sp_attrs:
        ast = ASTNode(AST_NAME)
        ast.setName(str(sp))
        all_species[ str(sp) ] = ("sp{}".format(sp_index), attr, ast)
        sp_index += 1
    #import ipdb; ipdb.set_trace()
    for sp in nw_model.list_species():
        if not all_species.has_key( sp.name() ):
            ast = ASTNode(AST_NAME)
            ast.setName( sp.name() )
            all_species[ sp.name() ] = ("sp{}".format(sp_index), 0, ast)
            sp_index += 1
    # Add to SBML model object
    for (sp, (sp_id, attr, ast) ) in all_species.items():
        sbml_sp = model.createSpecies()
        sbml_sp.setId(sp_id)
        sbml_sp.setName(str(sp))
        if isinstance(attr, (int, long, float, complex)):
            sbml_sp.setInitialAmount(attr)
        sbml_sp.setCompartment(compName)
    #}}}

    # Draw Reactions {{{
    astCytosol = ASTNode(AST_NAME)
    astCytosol.setName("cytosol");
    r_index = 0
    for rr in nw_model.reaction_rules():
        sbml_reaction = model.createReaction()
        sbml_reaction.setId("r{}".format(r_index))
        kl = sbml_reaction.createKineticLaw()
        multiple_factors = list()
        for reactant in rr.reactants():
            sbml_spr = sbml_reaction.createReactant()
            (sbml_spid, attr, ast) = all_species[ reactant.name()]
            #sbml_spr.setSpecies( all_species[ reactant.name() ][0] )
            sbml_spr.setSpecies( sbml_spid)
            multiple_factors.append(ast)
        for product in rr.products():
            sbml_spr = sbml_reaction.createProduct()
            (sbml_spid, attr, ast) = all_species[ product.name() ]
            sbml_spr.setSpecies( sbml_spid)
        # Kinetic Law related
        astKon = ASTNode(AST_NAME)
        astKon.setName("k{}".format(r_index))
        #multiple_factors.append(astKon)
        # Build Tree
        #import ipdb; ipdb.set_trace()
        print multiple_factors
        ast_number = 1
        '''
        while 1 < len(multiple_factors):
            astTimes_top = ASTNode(AST_TIMES)
            astTimes_top.addChild( multiple_factors.pop() )
            astTimes_top.addChild( multiple_factors.pop() )
            multiple_factors.insert(0, astTimes_top)
        kl.setMath( multiple_factors[0] )
        '''
        #asttimes = ASTNode(AST_TIMES)
        Node_top = ASTNode(AST_TIMES)
        kl.setMath(Node_top)
        Node_top.addChild(astKon)
        current_Node = Node_top

        import ipdb; ipdb.set_trace()
        while multiple_factors:
            if len(multiple_factors) == 1:
                current_Node.addChild(multiple_factors.pop(0))
            elif 1 < len(multiple_factors):
                child_Node = ASTNode(AST_TIMES)
                child_Node.addChild(multiple_factors.pop(0))
                current_Node.addChild(child_Node)
                current_Node = child_Node
            else:
                ast_one = ASTNode(AST_REAL) 
                current_Node.addChild(ast_one)
        r_index += 1
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

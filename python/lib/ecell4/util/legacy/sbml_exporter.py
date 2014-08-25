
import sys
import os.path
from libsbml import *
import ecell4.core as core
import ecell4.util.legacy.species 

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

        # Kinetic Parameter
        astKon = ASTNode(AST_NAME)
        astKon.setName("k{}".format(r_index))

        multiple_factors = []
        multiple_factors.append( astCytosol.deepCopy() )
        multiple_factors.append(astKon)
        for reactant in rr.reactants():
            sbml_spr = sbml_reaction.createReactant()
            (sbml_spid, attr, ast) = all_species[ reactant.name()]
            #sbml_spr.setSpecies( all_species[ reactant.name() ][0] )
            sbml_spr.setSpecies( sbml_spid)
            multiple_factors.append( ast.deepCopy() )
        for product in rr.products():
            sbml_spr = sbml_reaction.createProduct()
            (sbml_spid, attr, ast) = all_species[ product.name() ]
            sbml_spr.setSpecies( sbml_spid)
        # Build Tree
        #asttimes = ASTNode(AST_TIMES)
        '''
        if 2 < len(multiple_factors):
            current_node = ASTNode(AST_TIMES)
            current_node.addChild( multiple_factors.pop(0))
            kl.setMath(current_node)
            while 0 < len(multiple_factors):
                if len(multiple_factors) == 1:
                    current_node.addChild(multiple_factors.pop(0))
                elif 1 < len(multiple_factors):
                    new_node = ASTNode(AST_TIMES)
                    new_node.addChild(multiple_factors.pop(0))
                    current_node.addChild(new_node)
                    current_node = new_node
        '''
        import ipdb;ipdb.set_trace()
        ast2 = ASTNode(AST_TIMES)
        ast2.addChild(astKon)
        if len( rr.reactants() ) == 1:
            # 1 molecule reaction
            (sbml_spid, attr, ast_reactant) = all_species[rr.reactants()[0].name()]
            ast2.addChild( ast_reactant.deepCopy())

        elif len( rr.reactants() ) == 2:
            (sbml_spid0, attr0, ast_reactant0) = all_species[rr.reactants()[0].name()]
            (sbml_spid1, attr1, ast_reactant1) = all_species[rr.reactants()[1].name()]
            # 2 molecule reaction
            ast3 = ASTNode(AST_TIMES)
            ast3.addChild( ast_reactant0.deepCopy() )
            ast3.addChild( ast_reactant1.deepCopy() )
            ast2.addChild(ast3)
        else:
            raise RuntimeError("substances are too many.")
        ast1 = ASTNode(AST_TIMES)
        ast1.addChild(astCytosol)
        ast1.addChild(ast2)
        kl.setMath(ast1)
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

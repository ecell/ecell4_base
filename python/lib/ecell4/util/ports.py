import itertools
import ecell4


def export_sbml(model, y0={}, volume=1.0):
    """
    Export a model as a SBMLDocument.

    Parameters
    ----------
    model : NetworkModel or ODENetworkModel
    y0 : dict
        Initial condition.
    volume : Real or Real3, optional
        A size of the simulation volume.

    """
    import libsbml

    document = libsbml.SBMLDocument(3, 1)
    ns = libsbml.XMLNamespaces()
    ns.add("http://www.ecell.org/ns/ecell4", "ecell4")  #XXX: DUMMY URI
    document.setNamespaces(ns)
    m = document.createModel()

    comp1 = m.createCompartment()
    comp1.setId('world')
    comp1.setConstant(True)
    if isinstance(volume, ecell4.Real3):
        comp1.setSize(volume[0] * volume[1] * volume[2])
    else:
        comp1.setSize(volume)
    comp1.setSpatialDimensions(3)

    species_list = []
    for rr in model.reaction_rules():
        for sp in itertools.chain(rr.reactants(), rr.products()):
            species_list.append(sp)
    species_list = list(set(species_list))
    species_list.sort()

    for sp in species_list:
        s1 = m.createSpecies()
        s1.setId(sp.serial())
        s1.setCompartment('world')
        s1.setConstant(False)
        if sp.serial() in y0.keys():
            s1.setInitialAmount(y0[sp.serial()])
        else:
            s1.setInitialAmount(0)

        s1.setBoundaryCondition(False)
        s1.setHasOnlySubstanceUnits(False)

        s1.appendAnnotation('<annotation><ecell4:extension><ecell4:species serial="{:s}"/></ecell4:extension></annotation>'.format(sp.serial()))

    if isinstance(model, (ecell4.NetworkModel, ecell4.Model)):
        for cnt, rr in enumerate(model.reaction_rules()):
            r1 = m.createReaction()
            r1.setId("r{:d}".format(cnt))
            r1.setReversible(False)
            r1.setFast(False)

            kinetic_law = r1.createKineticLaw()
            # p1 = kinetic_law.createLocalParameter()
            # p1.setId("k")
            p1 = m.createParameter()
            p1.setId("k{:d}".format(cnt))
            p1.setConstant(True)
            p1.setValue(rr.k())

            species_coef_map = {}
            for sp in rr.reactants():
                if sp not in species_coef_map.keys():
                    species_coef_map[sp] = 1
                else:
                    species_coef_map[sp] += 1

            # math_exp = "k"
            math_exp = "k{:d}".format(cnt)
            for sp, coef in species_coef_map.items():
                s1 = r1.createReactant()
                s1.setSpecies(sp.serial())
                s1.setConstant(False)
                s1.setStoichiometry(coef)
                if coef == 1:
                    math_exp += "*{:s}".format(sp.serial())
                else:
                    math_exp += "*pow({:s},{:g})".format(sp.serial(), coef)

            species_coef_map = {}
            for sp in rr.products():
                if sp not in species_coef_map.keys():
                    species_coef_map[sp] = 1
                else:
                    species_coef_map[sp] += 1

            for sp, coef in species_coef_map.items():
                s1 = r1.createProduct()
                s1.setSpecies(sp.serial())
                s1.setConstant(False)
                s1.setStoichiometry(coef)

            math_ast = libsbml.parseL3Formula(math_exp)
            kinetic_law.setMath(math_ast)

    elif isinstance(model, ecell4.ode.ODENetworkModel):
        for cnt, rr in enumerate(model.reaction_rules()):
            r1 = m.createReaction()
            r1.setId("r{:d}".format(cnt))
            r1.setReversible(True)
            r1.setFast(False)

            kinetic_law = r1.createKineticLaw()

            species_coef_map = {}
            for sp, coef in zip(rr.reactants(), rr.reactants_coefficients()):
                if sp not in species_coef_map.keys():
                    species_coef_map[sp] = coef
                else:
                    species_coef_map[sp] += coef

            if rr.is_massaction():
                p1 = m.createParameter()
                p1.setId("k{:d}".format(cnt))
                # p1 = kinetic_law.createLocalParameter()
                # p1.setId("k")
                p1.setConstant(True)
                p1.setValue(rr.k())
                # math_exp = "k"
                math_exp = "k{:d}".format(cnt)
                for sp, coef in species_coef_map.items():
                    if coef == 1.0:
                        math_exp += "*{:s}".format(sp.serial())
                    else:
                        math_exp += "*pow({:s},{:g})".format(sp.serial(), coef)
            else:
                math_exp = rr.get_ratelaw().as_string()

            for sp, coef in species_coef_map.items():
                s1 = r1.createReactant()
                s1.setSpecies(sp.serial())
                s1.setConstant(False)
                s1.setStoichiometry(coef)

            species_coef_map = {}
            for sp in rr.products():
                if sp not in species_coef_map.keys():
                    species_coef_map[sp] = 1
                else:
                    species_coef_map[sp] += 1

            for sp, coef in species_coef_map.items():
                s1 = r1.createProduct()
                s1.setSpecies(sp.serial())
                s1.setConstant(False)
                s1.setStoichiometry(coef)

            math_ast = libsbml.parseL3Formula(math_exp)
            kinetic_law.setMath(math_ast)

    document.validateSBML()
    num_errors = (document.getNumErrors(libsbml.LIBSBML_SEV_ERROR)
                  + document.getNumErrors(libsbml.LIBSBML_SEV_FATAL))
    if num_errors > 0:
        messages = "The generated document is not valid."
        messages += " {} errors were found:\n".format(num_errors)
        for i in range(document.getNumErrors(libsbml.LIBSBML_SEV_ERROR)):
            err = document.getErrorWithSeverity(i, libsbml.LIBSBML_SEV_ERROR)
            messages += "{}: {}\n".format(err.getSeverityAsString(), err.getShortMessage())
        for i in range(document.getNumErrors(libsbml.LIBSBML_SEV_FATAL)):
            err = document.getErrorWithSeverity(i, libsbml.LIBSBML_SEV_FATAL)
            messages += "{}: {}\n".format(err.getSeverityAsString(), err.getShortMessage())
        raise RuntimeError(messages)

    return document

def save_sbml(filename, model, y0={}, volume=1.0):
    """
    Save a model in the SBML format.

    Parameters
    ----------
    model : NetworkModel or ODENetworkModel
    y0 : dict
        Initial condition.
    volume : Real or Real3, optional
        A size of the simulation volume.

    """
    import libsbml

    document = export_sbml(model, y0, volume)

    # with open(filename, 'w') as fout:
    #     fout.write(libsbml.writeSBMLToString(document))
    # writer = libsbml.SBMLWriter()
    # writer.writeSBML(document, filename)
    libsbml.writeSBML(document, filename)

    # reader = libsbml.SBMLReader()
    # document = reader.readSBML(filename)
    # if document.getNumErrors() > 0:
    #     document.printErrors()

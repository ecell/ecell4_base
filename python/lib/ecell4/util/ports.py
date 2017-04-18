import itertools
import copy
import warnings
import re

import ecell4

def replace_parseobj(expr, substitutes=None):
    substitutes = substitutes or {}

    import ecell4.util.decorator_base
    obj = ecell4.util.decorator_base.just_parse().evaluate(expr)

    from ecell4.util.decorator import traverse_ParseObj
    keys = []
    newexpr = str(traverse_ParseObj(copy.deepcopy(obj), keys))
    names = []
    for key in keys:
        if key in substitutes.keys():
            names.append(substitutes[key])
        else:
            raise RuntimeError(
                'unknown variable [{}] was used.'.format(key))
    return newexpr.format(*names)

def export_sbml(model, y0=None, volume=1.0):
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
    y0 = y0 or {}

    import libsbml

    document = libsbml.SBMLDocument(3, 1)

    # ns = libsbml.XMLNamespaces()
    # ns.add("http://www.ecell.org/ns/ecell4", "ecell4")  #XXX: DUMMY URI
    # document.setNamespaces(ns)

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

    sid_map = {}
    for cnt, sp in enumerate(species_list):
        sid_map[sp.serial()] = "s{:d}".format(cnt)

    for sp in species_list:
        sid = sid_map[sp.serial()]
        s1 = m.createSpecies()
        s1.setId(sid)
        s1.setName(sp.serial())
        s1.setCompartment('world')
        s1.setConstant(False)
        if sp.serial() in y0.keys():
            s1.setInitialAmount(y0[sp.serial()])
        else:
            s1.setInitialAmount(0)

        s1.setBoundaryCondition(False)
        s1.setHasOnlySubstanceUnits(False)

        # s1.appendAnnotation('<annotation><ecell4:extension><ecell4:species serial="{:s}"/></ecell4:extension></annotation>'.format(sp.serial()))

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
                sid = sid_map[sp.serial()]
                s1 = r1.createReactant()
                s1.setSpecies(sid)
                s1.setConstant(False)
                s1.setStoichiometry(coef)
                if coef == 1:
                    math_exp += "*{:s}".format(sid)
                else:
                    math_exp += "*pow({:s},{:g})".format(sid, coef)

            species_coef_map = {}
            for sp in rr.products():
                if sp not in species_coef_map.keys():
                    species_coef_map[sp] = 1
                else:
                    species_coef_map[sp] += 1

            for sp, coef in species_coef_map.items():
                sid = sid_map[sp.serial()]
                s1 = r1.createProduct()
                s1.setSpecies(sid)
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
                    sid = sid_map[sp.serial()]
                    if coef == 1.0:
                        math_exp += "*{:s}".format(sid)
                    else:
                        math_exp += "*pow({:s},{:g})".format(sid, coef)
            else:
                math_exp = rr.get_ratelaw().as_string()
                if math_exp in ('', '<lambda>'):
                    warnings.warn(
                        "The given ODEReactionRule [{:s}] might be invalid.".format(
                            rr.as_string()))
                math_exp = replace_parseobj(math_exp, sid_map)

            for sp, coef in species_coef_map.items():
                sid = sid_map[sp.serial()]
                s1 = r1.createReactant()
                s1.setSpecies(sid)
                s1.setConstant(False)
                s1.setStoichiometry(coef)

            species_coef_map = {}
            for sp, coef in zip(rr.products(), rr.products_coefficients()):
                if sp not in species_coef_map.keys():
                    species_coef_map[sp] = coef
                else:
                    species_coef_map[sp] += coef

            for sp, coef in species_coef_map.items():
                sid = sid_map[sp.serial()]
                s1 = r1.createProduct()
                s1.setSpecies(sid)
                s1.setConstant(False)
                s1.setStoichiometry(coef)

            math_ast = libsbml.parseL3Formula(math_exp)
            kinetic_law.setMath(math_ast)

    else:
        raise ValueError(
            "The invalid type of a Model was given [{:s}].".format(str(model))
            + " NetworkModel or ODENetworkModel must be given.")

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

def save_sbml(filename, model, y0=None, volume=1.0):
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
    y0 = y0 or {}

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

def import_sbml(document):
    """
    Import a model from a SBMLDocument.

    Parameters
    ----------
    document : SBMLDocument

    Returns
    -------
    model : NetworkModel or ODENetworkModel
    y0 : dict
        Initial condition.
    volume : Real or Real3, optional
        A size of the simulation volume.

    """
    from ecell4.util.decorator import generate_ratelaw

    m = document.getModel()

    if m.getNumCompartments() == 0:
        raise RuntimeError("No compartment was found.")
    elif m.getNumCompartments() > 1:
        warnings.warn(
            "[{:d}] compartments were found.".format(m.getNumCompartments())
            + " The second or later ones would be omitted.")

    comp1 = m.getCompartment(0)
    volume = comp1.getVolume()

    y0 = {}
    sid_map = {}
    for s1 in m.getListOfSpecies():
        sid = s1.getId()
        serial = s1.getName()
        sid_map[sid] = serial
        value = s1.getInitialAmount()
        if value != 0:
            y0[serial] = value

    kmap = {}
    for p1 in m.getListOfParameters():
        pid = p1.getId()
        if not re.match("^k[0-9]+$", pid):
            warnings.warn(
                "Parameter [{:s}] was just ommited.".format(pid))
        rid = "r{:s}".format(pid[1: ])
        kmap[rid] = p1.getValue()

    is_ode = False
    rrs = []

    for r1 in m.getListOfReactions():
        rid = r1.getId()

        is_massaction = (rid in kmap.keys())
        if is_massaction:
            k = kmap[rid]
        else:
            kinetic_law = r1.getKineticLaw()
            formula = kinetic_law.getFormula()
            k = replace_parseobj(formula, sid_map)

        reactants, products = [], []

        #XXX: The order of reactants is not consistent
        for s1 in r1.getListOfReactants():
            sid = s1.getSpecies()
            if sid not in sid_map:
                raise RuntimeError(
                    "Unknown Species' Id [{:s}] was given".format(sid))
            serial = sid_map[sid]
            coef = s1.getStoichiometry()
            reactants.append((serial, coef))

        #XXX: The order of products is not consistent
        for s1 in r1.getListOfProducts():
            sid = s1.getSpecies()
            if sid not in sid_map:
                raise RuntimeError(
                    "Unknown Species' Id [{:s}] was given".format(sid))
            serial = sid_map[sid]
            coef = s1.getStoichiometry()
            products.append((serial, coef))

        if (not is_massaction
            or len(reactants) > 2
            or any([coef not in (1, 2) for sp, coef in reactants])
            or any([not coef.is_integer() for sp, coef in products])
            or (len(reactants) == 2 and (reactants[0][1] == 2 or reactants[1][1] == 2))):
            is_ode = True
            rr = ecell4.ode.ODEReactionRule()

            for serial, coef in reactants:
                rr.add_reactant(ecell4.Species(serial), coef)
            for serial, coef in products:
                rr.add_product(ecell4.Species(serial), coef)

            if is_massaction:
                rr.set_k(k)
            else:
                func = generate_ratelaw(k, rr)
                rr.set_ratelaw(ecell4.ode.ODERatelawCallback(func, k))
        else:
            if len(reactants) == 1 and reactants[0][1] == 2:
                reactants[0] = (reactants[0][0], 1)
                reactants.append(reactants[0])

            rr = ecell4.ReactionRule()
            for serial, coef in reactants:
                rr.add_reactant(ecell4.Species(serial))
            for serial, coef in products:
                for _ in range(int(coef)):
                    rr.add_product(ecell4.Species(serial))
            rr.set_k(k)

        rrs.append(rr)

    m = ecell4.ode.ODENetworkModel() if is_ode else ecell4.NetworkModel()
    for rr in rrs:
        m.add_reaction_rule(rr)

    return m, y0, volume

def load_sbml(filename):
    """
    Load a model from a SBML file.

    Parameters
    ----------
    filename : str
        The input SBML filename.

    Returns
    -------
    model : NetworkModel or ODENetworkModel
    y0 : dict
        Initial condition.
    volume : Real or Real3, optional
        A size of the simulation volume.

    """
    import libsbml

    document = libsbml.readSBML(filename)
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
    return import_sbml(document)

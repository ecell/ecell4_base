import libsbml

try:
    from urllib2 import Request, urlopen
except ImportError:
    from urllib.request import Request, urlopen

class LambdaFunction(object):

    def __init__(self, args, formula, evalfunc):
        self.args = args
        self.formula = formula
        self.evalfunc = evalfunc

    def __call__(self, *args):
        return self.evalfunc(self.formula, dict(zip(self.args, args)))

class SBMLDataSource(object):

    def __init__(self, filename=None):
        if filename is not None:
            self.read(filename)

    def function_definitions(self, evalfunc=None, kwargs=None):
        kwargs = kwargs or {}
        for func in self.model.function_definitions:
            args = [func.getArgument(i).getName() for i in range(func.getNumArguments())]
            formula = libsbml.formulaToString(func.getBody())
            if evalfunc is None:
                yield (func.id, (args, formula))
            else:
                yield (func.id, LambdaFunction(args, formula, evalfunc))

    def read(self, filename):
        self.data = libsbml.SBMLReader().readSBML(filename)
        self.model = self.data.getModel()

    def read_from_string(self, xml):
        try:
            if isinstance(xml, unicode):
                xml = xml.encode('utf-8')
        except NameError:
            pass  # Python3
        self.data = libsbml.SBMLReader().readSBMLFromString(xml)
        self.model = self.data.getModel()

    def initial_amounts(self):
        for sp in self.model.species:
            if sp.isSetInitialAmount():
                yield (sp.id, sp.initial_amount)

    def initial_concentration(self):
        for sp in self.model.species:
            if sp.isSetInitialConcentration():
                yield (sp.id, sp.initial_concentration)

    def compartments(self):
        for comp in self.model.compartments:
            yield (comp.id, comp.volume)

    def parameters(self):
        for p in self.model.parameters:
            yield (p.id, p.value)

    def constants(self):
        for sp in self.model.species:
            if sp.getConstant():
                yield (sp.id)

    def assignment_rules(self, evalfunc=None, kwargs=None):
        kwargs = kwargs or {}

        for rule in self.model.rules:
            if rule.isAssignment():
                if evalfunc is None:
                    yield (rule.variable, rule.formula)
                else:
                    #XXX: Why not evaluate variable?
                    yield (rule.variable, evalfunc(rule.formula, kwargs))

    def reactions(self, evalfunc=None, kwargs=None):
        kwargs = kwargs or {}

        for r in self.model.reactions:
            reactants = [(reactant.species, reactant.stoichiometry)
                         for reactant in r.reactants]
            products = [(product.species, product.stoichiometry)
                        for product in r.products]

            formula = r.getKineticLaw().formula
            params = dict((p.id, p.value) for p in r.getKineticLaw().parameters)

            if evalfunc is None:
                yield (reactants, products, formula, params)
            else:
                params.update(kwargs)
                yield (sum((evalfunc(sp, params) * coef for sp, coef in reactants), evalfunc('~EmptySet')),
                       sum((evalfunc(sp, params) * coef for sp, coef in products), evalfunc('~EmptySet')),
                       evalfunc(formula, params))

class BioModelsDataSource(SBMLDataSource):

    URL = "https://www.ebi.ac.uk/biomodels-main/download?mid={mid}"
    XML = {}

    def __init__(self, mid, cache=True):
        SBMLDataSource.__init__(self)

        if not cache or mid not in self.XML.keys():
            url = self.URL.format(mid=mid)
            req = Request(url)
            response = urlopen(req)
            xml = response.read().decode('utf-8')
            if cache:
                self.XML[mid] = xml
        else:
            xml = self.XML[mid]

        #import pdb; pdb.set_trace()
        self.read_from_string(xml)

# from bs4 import BeautifulSoup
# import requests
# 
# def biomodels(biomodels_id):
#     xml_soup = BeautifulSoup(requests.get("https://www.ebi.ac.uk/biomodels-main/download?mid=" + biomodels_id).content, 'xml')
#     reactions = xml_soup.listOfReactions.find_all('reaction')
#     rpks = []
#     for r in reactions:
#         rpk = []
#         reactant = r.listOfReactants.find('speciesReference')['species']
#         product = r.listOfProducts.find('speciesReference')['species']
#         kinetics = {}
#         for k in r.listOfParameters.find_all('parameter'):
#             kinetics[k['id']] = k['value']
#         rpks.append([reactant, product, kinetics])
#     return rpks


if __name__ == '__main__':
    import sys

    from ecell4 import *
    biomodels = BioModelsDataSource

    mid = 'BIOMD0000000005'

    y0 = dict(biomodels(mid).initial_amounts())
    print(y0)

    params = dict(biomodels(mid).parameters())
    params.update(biomodels(mid).compartments())
    print(params)


    with reaction_rules():
        params['EmptySet'] = ~EmptySet  #XXX: Just ignore EmptySet
        params.update(dict(biomodels(mid).assignment_rules(_eval, params)))
        params.update(dict(biomodels(mid).function_definitions(_eval, params)))

        for reactants, products, formula in biomodels(mid).reactions(_eval, params):
            reactants > products | formula

    m = get_model()
    print([rr.as_string() for rr in m.reaction_rules()])

    print(run_simulation(100, model=m, y0=y0, return_type='array'))

    # sbml = SBMLDataSource
    #     for reactants, products, formula, parameters in sbml(filename).reactions():

    # filename = sys.argv[1]

    # y0 = dict(sbml(filename).initial_amounts())
    # # y0.update(dict(sbml(filename).compartments()))
    # print(y0)

    # params = dict(sbml(filename).parameters())
    # params.update(dict(sbml(filename).compartments()))
    # params['compartment'] = 1.0
    # print(params)

    # with reaction_rules():
    #     # params.update(dict((var, _eval(formula)) for var, formula in sbml(filename).assignment_rules()))
    #     print(dict((var, _eval(formula)) for var, formula in sbml(filename).assignment_rules()))

    #     for reactants, products, formula, parameters in sbml(filename).reactions():
    #         parameters.update(params)

    #         (sum((_eval(sp) * coef for sp, coef in reactants), ~EmptySet)
    #                 > sum((_eval(sp) * coef for sp, coef in products), ~EmptySet) | _eval(formula, parameters))

    #         # _sum(reactants) > _sum(products) | _eval(formula, params)
    #         # _sum(_mul(reactants, reactant_coefficients)) > _sum(_mul(products, product_coefficients)) | _eval(formula, params)

    # m = get_model()
    # print([rr.as_string() for rr in m.reaction_rules()])

    # run_simulation(60, model=m, y0=y0, opt_kwargs={'legend': False})
    # # w = run_simulation(0.0032, model=m, y0=y0, species_list=['x1'], return_type='world')
    # # y0 = dict((sp.serial(), w.get_value(sp)) for sp in w.list_species())
    # # y0['k5'] = 1.55
    # # run_simulation(0.1, model=m, y0=y0, species_list=['x1'])

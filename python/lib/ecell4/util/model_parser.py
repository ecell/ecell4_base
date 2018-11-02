import copy
import numbers
import itertools
import math

# import operator
# from functools import reduce

from . import parseobj
# from ..extra import unit

import ecell4.core

ENABLE_RATELAW = True
ENABLE_IMPLICIT_DECLARATION = True

RATELAW_RESERVED_FUNCTIONS = {
    'pow': pow, 'exp': math.exp, 'log': math.log,
    'sin': math.sin, 'cos': math.cos, 'tan': math.tan,
    'asin': math.asin, 'acos': math.acos, 'atan': math.atan
    }

RATELAW_RESERVED_CONSTANTS = {
    '_t': None,  #XXX: just reserved
    'pi': math.pi
    }

def define_ratelaw_function(key, val):
    global RATELAW_RESERVED_FUNCTIONS
    if not callable(val):
        raise ValueError('A function must be callable [{}]'.format(type(val)))
    RATELAW_RESERVED_FUNCTIONS[key] = val

def define_ratelaw_constant(key, val):
    global RATELAW_RESERVED_CONSTANTS
    if not isinstance(val, numbers.Number):
        raise ValueError('A constant must be a number [{}]'.format(type(val)))
    RATELAW_RESERVED_CONSTANTS[key] = val

def generate_species(obj):
    if not isinstance(obj, (parseobj.AnyCallable, parseobj.ParseObj)):
        raise TypeError(
            "Argument 1 must be AnyCallable or ParseObj."
            " '{}' was given [{}].".format(type(obj).__name__, obj))
    return ecell4.core.Species(str(obj))

def generate_species_with_coefficient(obj):
    if isinstance(obj, (parseobj.AnyCallable, parseobj.ParseObj)):
        return (generate_species(obj), None)
    elif isinstance(obj, parseobj.InvExp):
        return (generate_species(obj._target), 0)
    elif isinstance(obj, parseobj.MulExp):
        sp, coef = None, 1.0
        for elem in obj._elements():
            if isinstance(elem, numbers.Real):
                coef *= elem
                continue
            if sp is not None:
                raise RuntimeError(
                    "Species is already assigned [{}].".format(sp.serila())
                    + " [{!r}] was given duplicately.".format(elem))
            sp = generate_species(elem)
        return (sp, coef)
    else:
        raise TypeError(
            "Argument 1 must be AnyCallable, ParseObj, InvExp or MulExp."
            " '{}' was given [{}].".format(type(obj).__name__, obj))

def generate_list_of_species_with_coefficient(obj):
    if isinstance(obj, parseobj.AddExp):
        return [generate_species_with_coefficient(elem) for elem in obj._elements()]
    else:
        return [generate_species_with_coefficient(obj)]

def generate_reaction_rule_options(elements):
    if elements is None or len(elements) == 0:
        return {}

    opts = {}
    for elem in elements:
        if isinstance(elem, ecell4.core.ReactionRulePolicy):
            if 'policy' not in opts.keys():
                opts['policy'] = elem.get()
            else:
                opts['policy'] |= elem.get()
        else:
            if 'k' in opts.keys():
                raise RuntimeError(
                    "A kinetic rate is already assigned [{}].".format(opts['k'])
                    + " [{}] was given duplicately.".format(elem))
            opts['k'] = elem
    return opts

def generate_reaction_rule(lhs, rhs, k=None, policy=None):
    if k is None:
        raise RuntimeError('A kinetic rate must be given.')

    rr = ecell4.core.ReactionRule([sp for (sp, _) in lhs], [sp for (sp, _) in rhs])

    if (callable(k)  # Function
            or (ENABLE_RATELAW and isinstance(k, (parseobj.ExpBase, parseobj.AnyCallable)))  # Formula
            or any([coef is not None for (_, coef) in itertools.chain(lhs, rhs)])):  # Stoichiometry

        if ENABLE_RATELAW and isinstance(k, (parseobj.ExpBase, parseobj.AnyCallable)):
            func, name = generate_ratelaw(k, rr, ENABLE_IMPLICIT_DECLARATION)
            desc = ecell4.core.ReactionRuleDescriptorPyfunc(func, name)
        elif callable(k):
            desc = ecell4.core.ReactionRuleDescriptorPyfunc(k, "")
        else:
            if not isinstance(k, (numbers.Real, ecell4.core.Quantity)):
                raise TypeError(
                    "A kinetic rate must be float or Quantity."
                    "'{}' was given [{}].".format(type(k).__name__, k))
            desc = ecell4.core.ReactionRuleDescriptorMassAction(k)

        desc.set_reactant_coefficients([coef or 1 for (_, coef) in lhs])
        desc.set_product_coefficients([coef or 1 for (_, coef) in rhs])

        rr.set_descriptor(desc)
    elif isinstance(k, (numbers.Real, ecell4.core.Quantity)):
        rr.set_k(k)
    # elif unit.HAS_PINT and isinstance(k, unit._Quantity):  # Kinetic rate given as a quantity
    #     if unit.STRICT:
    #         if len(lhs) == 0 and not unit.check_dimensionality(k, '1/[time]/[volume]'):
    #             raise ValueError(
    #                 "Cannot convert [k] from '{}' ({}) to '1/[time]/[volume]'".format(k.dimensionality, k.u))
    #         elif not unit.check_dimensionality(k, '1/[time]' + '/[concentration]' * (len(lhs) - 1)):
    #             raise ValueError(
    #                 "Cannot convert [k] from '{}' ({}) to '{}'".format(
    #                     k.dimensionality, k.u, '1/[time]' + '/[concentration]' * (len(lhs) - 1)))
    #     rr = ecell4.core.ReactionRule([sp for (sp, _) in lhs], [sp for (sp, _) in rhs], k.to_base_units().magnitude)
    #     if policy is not None:
    #         rr.set_policy(policy)
    #     return rr
    else:
        raise TypeError(
            "A kinetic rate must be float, Quantity or function."
            " '{}' was given [{}].".format(type(k).__name__, k))

    if policy is not None:
        if not isinstance(policy, ecell4.core.ReactionRulePolicy):
            raise TypeError(
                "policy must be ReactionRulePolicy."
                " '{}' was given [{}].".format(type(policy).__name__, policy))
        rr.set_policy(policy)

    return rr

class Visitor(object):

    def visit_const(self, obj):
        return self.visit_default(obj)

    def visit_species(self, obj):
        return self.visit_default(obj)

    def visit_quantity(self, obj):
        return self.visit_default(obj)

    def visit_func(self, obj, *args):
        return self.visit_default(obj)

    def visit_expression(self, obj, *args):
        return self.visit_default(obj)

    def visit_default(self, obj):
        return obj

def dispatch(obj, visitor):
    if isinstance(obj, parseobj.AnyCallable):
        obj = obj._as_ParseObj()

    if isinstance(obj, parseobj.ParseObj):
        if obj._size() == 1 and obj._elems[0].name in RATELAW_RESERVED_FUNCTIONS:
            # function
            subobj = obj._elems[0]
            assert subobj.key is None
            assert subobj.modification is None
            assert (subobj.args is not None and subobj.kwargs == {}) or (subobj.args is None and subobj.kwargs is None)
            args = [dispatch(subobj.args[i], visitor) for i in range(len(subobj.args))]
            return visitor.visit_func(obj, *args)
        elif obj._size() == 1 and obj._elems[0].name in RATELAW_RESERVED_CONSTANTS:
            # constant
            subobj = obj._elems[0]
            assert subobj.key is None
            assert subobj.modification is None
            assert subobj.args is None
            assert subobj.kwargs is None
            return visitor.visit_const(obj)
        else:
            # species
            return visitor.visit_species(obj)
    elif isinstance(obj, parseobj.ExpBase):
        args = [dispatch(obj._elems[i], visitor) for i in range(len(obj._elems))]
        return visitor.visit_expression(obj, *args)
    # elif isinstance(obj, unit._Quantity):
    #     return visitor.visit_quantity(obj)
    else:
        return visitor.visit_default(obj)

class SpeciesParsingVisitor(Visitor):

    def __init__(self):
        Visitor.__init__(self)
        self.__keys = []
        # self.__quantities = []

    @property
    def keys(self):
        return self.__keys

    # @property
    # def quantities(self):
    #     return self.__quantities

    def visit_species(self, obj):
        serial = ecell4.core.Species(str(obj)).serial()
        if serial in self.__keys:
            return "{{{0:d}}}".format(self.__keys.index(serial))
        self.__keys.append(serial)
        return "{{{0:d}}}".format(len(self.__keys) - 1)

    # def visit_quantity(self, obj):
    #     assert not isinstance(obj.magnitude, (parseobj.AnyCallable, parseobj.ExpBase))
    #     self.__quantities.append(obj)
    #     return obj.to_base_units().magnitude

    def visit_func(self, obj, *args):
        subobj = obj._elems[0]
        subobj.args = tuple(args)
        return obj

    def visit_expression(self, obj, *args):
        assert len(obj._elems) == len(args)
        obj._elems = list(args)
        return obj

# class DimensionalityCheckingVisitor(Visitor):
# 
#     OPERATORS = {
#         parseobj.PosExp: operator.pos,
#         parseobj.NegExp: operator.neg,
#         parseobj.SubExp: lambda *args: operator.sub, # reduce(operator.sub, args[1: ], args[0]),
#         parseobj.DivExp: operator.truediv, # lambda *args: reduce(operator.truediv, args[1: ], args[0]),
#         parseobj.PowExp: operator.pow,
#         parseobj.AddExp: lambda *args: reduce(operator.add, args[1: ], args[0]), # operator.add,
#         parseobj.MulExp: lambda *args: reduce(operator.mul, args[1: ], args[0]), # operator.mul,
#         # parseobj.InvExp: operator.inv,
#         # parseobj.AndExp: operator.and_,
#         # parseobj.GtExp: operator.gt,
#         # parseobj.NeExp: operator.ne,
#         # parseobj.EqExp: operator.eq,
#         }
# 
#     def __init__(self, ureg):
#         Visitor.__init__(self)
#         self.__ureg = ureg
# 
#     def visit_const(self, obj):
#         key = obj._elems[0].name
#         if key == '_t':
#             dim = self.__ureg.Quantity(1.0, "second").to_base_units().u
#         else:
#             dim = RATELAW_RESERVED_CONSTANTS[key]
#         return (obj, dim)
# 
#     def visit_species(self, obj):
#         dim = self.__ureg.Quantity(1.0, "molar").to_base_units().u
#         return (obj, dim)
# 
#     def visit_quantity(self, obj):
#         assert not isinstance(obj.magnitude, (parseobj.AnyCallable, parseobj.ExpBase))
#         val = obj.to_base_units()
#         val, dim = val.magnitude, val.u
#         return (val, dim)
# 
#     def visit_func(self, obj, *args):
#         func = RATELAW_RESERVED_FUNCTIONS[obj._elems[0].name]
#         val = func(*[1.0 * x for _, x in args])
#         dim = val.to_base_units().u
#         subobj = obj._elems[0]
#         subobj.args = tuple([x for x, _ in args])
#         return (obj, dim)
# 
#     def visit_expression(self, obj, *args):
#         assert len(obj._elems) == len(args)
#         for cls, op in self.OPERATORS.items():
#             if isinstance(obj, cls):
#                 val = op(*[1.0 * x for _, x in args])
#                 dim = val.to_base_units().u
#                 obj._elems = list([x for x, _ in args])
#                 return (obj, dim)
#         raise ValueError('Unknown dimensionality for the given object [{}]'.format(str(obj)))
# 
#     def visit_default(self, obj):
#         return (obj, obj)

def generate_ratelaw(obj, rr, implicit=False):
    label = str(obj)
    visitor = SpeciesParsingVisitor()
    exp = str(dispatch(copy.deepcopy(obj), visitor))

    # if unit.STRICT and len(visitor.quantities) > 0:
    #     ureg = visitor.quantities[0]._REGISTRY
    #     if any([q._REGISTRY != ureg for q in visitor.quantities[1: ]]):
    #         raise ValueError('Cannot operate with Quantity and Quantity of different registries.')
    #     label, ret = dispatch(copy.deepcopy(obj), DimensionalityCheckingVisitor(ureg))
    #     label = str(label)

    #     if not isinstance(ret, unit._Unit):
    #         ret = ureg.Unit('dimensionless')
    #     if not unit.check_dimensionality(ret, '[concentration]/[time]'):
    #         raise RuntimeError(
    #             "A rate law must have dimension '{}'. '{}' was given.".format(
    #                 ureg.get_dimensionality("[concentration]/[time]"), ret.dimensionality))

    aliases = {}
    for i, sp in enumerate(rr.reactants()):
        aliases[sp.serial()] = "_r[{0:d}]".format(i)
    for i, sp in enumerate(rr.products()):
        aliases[sp.serial()] = "_p[{0:d}]".format(i)
    names = []
    for key in visitor.keys:
        if key in aliases.keys():
            names.append(aliases[key])
        elif implicit:
            names.append("_r[{0:d}]".format(len(rr.reactants())))
            aliases[key] = names[-1]
            rr.add_reactant(ecell4.core.Species(key))
            rr.add_product(ecell4.core.Species(key))
        else:
            raise RuntimeError(
                'unknown variable [{}] was used.'.format(key))
    exp = exp.format(*names)
    # print(exp)
    f = eval("lambda _r, _p, _v, _t, _rc, _pc: {0}".format(exp))
    f.__globals__.update(RATELAW_RESERVED_FUNCTIONS)
    f.__globals__.update((key, val) for key, val in RATELAW_RESERVED_CONSTANTS if val is not None)

    return (f, label)
    # return (lambda _r, _p, *args: eval(exp))

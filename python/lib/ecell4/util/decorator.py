import copy
import types
import numbers
import warnings
import functools
import itertools
import operator
import math

from functools import reduce

from . import parseobj
from .decorator_base import Callback, JustParseCallback, ParseDecorator
# from ..extra import unit

import ecell4.core

ENABLE_RATELAW = True
ENABLE_IMPLICIT_DECLARATION = True

SPECIES_ATTRIBUTES = []
REACTION_RULES = []

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

def parse_species(obj):
    if not isinstance(obj, (parseobj.AnyCallable, parseobj.ParseObj)):
        raise TypeError(
            "Argument 1 must be AnyCallable or ParseObj."
            " '{}' was given [{}].".format(type(obj).__name__, obj))
    return ecell4.core.Species(str(obj))

class SpeciesAttributesCallback(Callback):

    def __init__(self):
        Callback.__init__(self)

        self.bitwise_operations = []

    def get(self):
        return copy.copy(self.bitwise_operations)

    def set(self):
        global SPECIES_ATTRIBUTES
        SPECIES_ATTRIBUTES.extend(self.bitwise_operations)

    def notify_bitwise_operations(self, obj):
        # attribute_dimensionality = {'D': '[length]**2/[time]', 'radius': '[length]'}

        if not isinstance(obj, parseobj.OrExp):
            raise TypeError('An invalid object was given [{}]'.format(repr(obj)))

        elems = obj._elements()
        rhs = elems[-1]
        if isinstance(rhs, parseobj.ExpBase):
            return
        elif not isinstance(rhs, dict):
            raise TypeError('parameter must be given as a dict; "{}" given'.format(type(rhs)))

        for lhs in elems[: -1]:
            sp = parse_species(lhs)

            for key, value in rhs.items():
                if not isinstance(key, str):
                    raise TypeError(
                        "Attribute key must be string."
                        " '{}' was given [{}].".format(type(key).__name__, key))
                if not isinstance(value, (numbers.Real, str, bool, ecell4.core.Quantity)):
                    raise TypeError(
                        "Attribute value must be int, float, string, boolean or Quantity."
                        " '{}' was given [{}].".format(type(value).__name__, value))

                sp.set_attribute(key, value)

                # if unit.HAS_PINT and isinstance(value, unit._Quantity):
                #     if (unit.STRICT and key in attribute_dimensionality
                #         and not unit.check_dimensionality(value, attribute_dimensionality[key])):
                #             raise ValueError("Cannot convert [{}] from '{}' ({}) to '{}'".format(
                #                 key, value.dimensionality, value.u, attribute_dimensionality[key]))
                #     sp.set_attribute(key, value.to_base_units().magnitude)
                # else:
                #     sp.set_attribute(key, value)

            self.bitwise_operations.append(sp)

    def notify_comparisons(self, obj):
        raise RuntimeError(
            'ReactionRule definitions are not allowed'
            + ' in "species_attributes"')

def parse_species_with_coefficient(obj):
    if isinstance(obj, (parseobj.AnyCallable, parseobj.ParseObj)):
        return (parse_species(obj), None)
    elif isinstance(obj, parseobj.InvExp):
        return (parse_species(obj._target), 0)
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
            sp = parse_species(elem)
        return (sp, coef)
    else:
        raise TypeError(
            "Argument 1 must be AnyCallable, ParseObj, InvExp or MulExp."
            " '{}' was given [{}].".format(type(obj).__name__, obj))

def parse_list_of_species_with_coefficient(obj):
    if isinstance(obj, parseobj.AddExp):
        return [parse_species_with_coefficient(elem) for elem in obj._elements()]
    else:
        return [parse_species_with_coefficient(obj)]

def parse_reaction_rule_options(elements):
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

def parse_reaction_rule(lhs, rhs, k=None, policy=None):
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

class ReactionRulesCallback(Callback):

    def __init__(self):
        Callback.__init__(self)

        self.comparisons = []

    def get(self):
        return copy.copy(self.comparisons)

    def set(self):
        global REACTION_RULES
        REACTION_RULES.extend(self.comparisons)

    def notify_comparisons(self, obj):
        if not isinstance(obj, (parseobj.EqExp, parseobj.GtExp, parseobj.LtExp)):
            raise RuntimeError('An ivalid operation was given [{!r}].'.format(obj))

        lhs, rhs = obj._lhs, obj._rhs

        if isinstance(obj._lhs, parseobj.OrExp):
            lhs = obj._lhs._elements()[0]
        else:
            lhs = obj._lhs

        if not isinstance(obj._rhs, parseobj.OrExp):
            raise RuntimeError(
                "A right-hand-side is ill-formed. OrExp must be given [{!r}].".format(obj._rhs))
        elif len(obj._rhs._elements()) == 0:
            raise RuntimeError('No product was given in the right-hand-side.')
        else:
            rhs = obj._rhs._elements()[0]
            opts = obj._rhs._elements()[1: ]

        reactants = parse_list_of_species_with_coefficient(lhs)
        products = parse_list_of_species_with_coefficient(rhs)
        reactants = tuple((sp, coef) for sp, coef in reactants if coef != 0)
        products = tuple((sp, coef) for sp, coef in products if coef != 0)

        opts = parse_reaction_rule_options(opts)
        if 'k' not in opts.keys():
            raise RuntimeError('No kinetic rate or law was given.')
        params = opts['k']

        if isinstance(obj, parseobj.EqExp):
            if not isinstance(params, (tuple, list)):
                raise RuntimeError(
                    "Parameter must be list or tuple."
                    " '{}' was given [{}].".format(type(params).__name__, params))
            elif len(params) != 2:
                raise RuntimeError(
                    "Parameter must have size, 2."
                    " '{}' was given [{}].".format(len(params), params))

            self.comparisons.append(
                parse_reaction_rule(reactants, products, params[0], opts.get('policy')))
            self.comparisons.append(
                parse_reaction_rule(products, reactants, params[1], opts.get('policy')))
        elif isinstance(obj, parseobj.GtExp):
            self.comparisons.append(
                parse_reaction_rule(reactants, products, params, opts.get('policy')))
        elif isinstance(obj, parseobj.LtExp):
            self.comparisons.append(
                parse_reaction_rule(products, reactants, params, opts.get('policy')))

def get_model(is_netfree=False, without_reset=False, seeds=None, effective=False):
    """
    Generate a model with parameters in the global scope, ``SPECIES_ATTRIBUTES``
    and ``REACTIONRULES``.

    Parameters
    ----------
    is_netfree : bool, optional
        Return ``NetfreeModel`` if True, and ``NetworkModel`` if else.
        Default is False.
    without_reset : bool, optional
        Do not reset the global variables after the generation if True.
        Default is False.
    seeds : list, optional
        A list of seed ``Species`` for expanding the model.
        If this is not None, generate a ``NetfreeModel`` once, and return a
        ``NetworkModel``, which is an expanded form of that with the given seeds.
        Default is None.
    effective : bool, optional
        See ``NetfreeModel.effective`` and ``Netfree.set_effective``.
        Only meaningfull with option ``is_netfree=True``.
        Default is False

    Returns
    -------
    model : NetworkModel, NetfreeModel

    """
    try:
        if seeds is not None or is_netfree:
            m = ecell4.core.NetfreeModel()
        else:
            m = ecell4.core.NetworkModel()

        for sp in SPECIES_ATTRIBUTES:
            m.add_species_attribute(sp)
        for rr in REACTION_RULES:
            m.add_reaction_rule(rr)

        if not without_reset:
            reset_model()

        if seeds is not None:
            return m.expand(seeds)

        if isinstance(m, ecell4.core.NetfreeModel):
            m.set_effective(effective)
    except Exception as e:
        reset_model()
        raise e

    return m

def reset_model():
    """
    Reset all values, ``SPECIES_ATTRIBUTES`` and ``REACTIONRULES``,
    in the global scope.

    """
    global SPECIES_ATTRIBUTES
    global REACTION_RULES

    SPECIES_ATTRIBUTES = []
    REACTION_RULES = []

reaction_rules = functools.partial(ParseDecorator, ReactionRulesCallback)
species_attributes = functools.partial(ParseDecorator, SpeciesAttributesCallback)

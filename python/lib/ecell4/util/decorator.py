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
from ..extra import unit

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

def generate_Species(obj):
    if isinstance(obj, parseobj.AnyCallable):
        obj = obj._as_ParseObj()

    if isinstance(obj, parseobj.ParseObj):
        return ((ecell4.core.Species(str(obj)), None), )
    elif isinstance(obj, parseobj.InvExp):
        return (None, )
    elif isinstance(obj, parseobj.MulExp):
        subobjs = obj._elements()

        retval, coef = None, 1
        for subobj in subobjs:
            if isinstance(subobj, numbers.Number):
                coef *= subobj
            elif retval is not None:
                raise RuntimeError(
                    'only a single species must be given; %s given'
                    % (repr(obj)))
            else:
                retval = generate_Species(subobj)

        return [(sp[0], coef if sp[1] is None else sp[1] * coef)
                if sp is not None else None
                    for sp in retval]

    elif isinstance(obj, parseobj.AddExp):
        subobjs = obj._elements()
        return tuple(itertools.chain(*[
            generate_Species(subobj) for subobj in subobjs]))
    else:
        raise RuntimeError('invalid expression; "%s" given' % str(obj))

def generate_ReactionRule(lhs, rhs, k=None):
    if k is None:
        raise RuntimeError('no parameter is specified')

    elif (callable(k)  # Function
          or (ENABLE_RATELAW
              and isinstance(k, (parseobj.ExpBase, parseobj.AnyCallable)))  # Formula
          or any([coef is not None
                  for (sp, coef) in itertools.chain(lhs, rhs)])):  # Stoichiometry
        from ecell4.ode import ODEReactionRule, ODERatelawCallback

        rr = ODEReactionRule()
        for sp, coef in lhs:
            rr.add_reactant(sp, coef or 1)
        for sp, coef in rhs:
            rr.add_product(sp, coef or 1)

        if ENABLE_RATELAW and isinstance(k, (parseobj.ExpBase, parseobj.AnyCallable)):
            func, name = generate_ratelaw(k, rr, ENABLE_IMPLICIT_DECLARATION)
            rr.set_ratelaw(ODERatelawCallback(func, name))
        elif callable(k):
            rr.set_ratelaw(ODERatelawCallback(k))
        else:
            rr.set_k(k)
        return rr

    elif isinstance(k, numbers.Number):  # Kinetic rate
        return ecell4.core.ReactionRule([sp for (sp, _) in lhs], [sp for (sp, _) in rhs], k)

    elif unit.HAS_PINT and isinstance(k, unit._Quantity):  # Kinetic rate given as a quantity
        if unit.STRICT:
            if len(lhs) == 0 and not unit.check_dimensionality(k, '1/[time]/[volume]'):
                raise ValueError(
                    "Cannot convert [k] from '{}' ({}) to '1/[time]/[volume]'".format(k.dimensionality, k.u))
            elif not unit.check_dimensionality(k, '1/[time]' + '/[concentration]' * (len(lhs) - 1)):
                raise ValueError(
                    "Cannot convert [k] from '{}' ({}) to '{}'".format(
                        k.dimensionality, k.u, '1/[time]' + '/[concentration]' * (len(lhs) - 1)))
        return ecell4.core.ReactionRule([sp for (sp, _) in lhs], [sp for (sp, _) in rhs], k.to_base_units().magnitude)

    raise ValueError('parameter must be given as a number; "%s" given' % str(k))

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
    elif isinstance(obj, unit._Quantity):
        return visitor.visit_quantity(obj)
    else:
        return visitor.visit_default(obj)

class SpeciesParsingVisitor(Visitor):

    def __init__(self):
        Visitor.__init__(self)
        self.__keys = []
        self.__quantities = []

    @property
    def keys(self):
        return self.__keys

    @property
    def quantities(self):
        return self.__quantities

    def visit_species(self, obj):
        serial = ecell4.core.Species(str(obj)).serial()
        if serial in self.__keys:
            return "{{{0:d}}}".format(self.__keys.index(serial))
        self.__keys.append(serial)
        return "{{{0:d}}}".format(len(self.__keys) - 1)

    def visit_quantity(self, obj):
        assert not isinstance(obj.magnitude, (parseobj.AnyCallable, parseobj.ExpBase))
        self.__quantities.append(obj)
        return obj.to_base_units().magnitude

    def visit_func(self, obj, *args):
        subobj = obj._elems[0]
        subobj.args = tuple(args)
        return obj

    def visit_expression(self, obj, *args):
        assert len(obj._elems) == len(args)
        obj._elems = list(args)
        return obj

class DimensionalityCheckingVisitor(Visitor):

    OPERATORS = {
        parseobj.PosExp: operator.pos,
        parseobj.NegExp: operator.neg,
        parseobj.SubExp: lambda *args: operator.sub, # reduce(operator.sub, args[1: ], args[0]),
        parseobj.DivExp: operator.truediv, # lambda *args: reduce(operator.truediv, args[1: ], args[0]),
        parseobj.PowExp: operator.pow,
        parseobj.AddExp: lambda *args: reduce(operator.add, args[1: ], args[0]), # operator.add,
        parseobj.MulExp: lambda *args: reduce(operator.mul, args[1: ], args[0]), # operator.mul,
        # parseobj.InvExp: operator.inv,
        # parseobj.AndExp: operator.and_,
        # parseobj.GtExp: operator.gt,
        # parseobj.NeExp: operator.ne,
        # parseobj.EqExp: operator.eq,
        }

    def __init__(self, ureg):
        Visitor.__init__(self)
        self.__ureg = ureg

    def visit_const(self, obj):
        key = obj._elems[0].name
        if key == '_t':
            dim = self.__ureg.Quantity(1.0, "second").to_base_units().u
        else:
            dim = RATELAW_RESERVED_CONSTANTS[key]
        return (obj, dim)

    def visit_species(self, obj):
        dim = self.__ureg.Quantity(1.0, "molar").to_base_units().u
        return (obj, dim)

    def visit_quantity(self, obj):
        assert not isinstance(obj.magnitude, (parseobj.AnyCallable, parseobj.ExpBase))
        val = obj.to_base_units()
        val, dim = val.magnitude, val.u
        return (val, dim)

    def visit_func(self, obj, *args):
        func = RATELAW_RESERVED_FUNCTIONS[obj._elems[0].name]
        val = func(*[1.0 * x for _, x in args])
        dim = val.to_base_units().u
        subobj = obj._elems[0]
        subobj.args = tuple([x for x, _ in args])
        return (obj, dim)

    def visit_expression(self, obj, *args):
        assert len(obj._elems) == len(args)
        for cls, op in self.OPERATORS.items():
            if isinstance(obj, cls):
                val = op(*[1.0 * x for _, x in args])
                dim = val.to_base_units().u
                obj._elems = list([x for x, _ in args])
                return (obj, dim)
        raise ValueError('Unknown dimensionality for the given object [{}]'.format(str(obj)))

    def visit_default(self, obj):
        return (obj, obj)

def generate_ratelaw(obj, rr, implicit=False):
    label = str(obj)
    visitor = SpeciesParsingVisitor()
    exp = str(dispatch(copy.deepcopy(obj), visitor))

    if unit.STRICT and len(visitor.quantities) > 0:
        ureg = visitor.quantities[0]._REGISTRY
        if any([q._REGISTRY != ureg for q in visitor.quantities[1: ]]):
            raise ValueError('Cannot operate with Quantity and Quantity of different registries.')
        label, ret = dispatch(copy.deepcopy(obj), DimensionalityCheckingVisitor(ureg))
        label = str(label)

        if not isinstance(ret, unit._Unit):
            ret = ureg.Unit('dimensionless')
        if not unit.check_dimensionality(ret, '[concentration]/[time]'):
            raise RuntimeError(
                "A rate law must have dimension '{}'. '{}' was given.".format(
                    ureg.get_dimensionality("[concentration]/[time]"), ret.dimensionality))

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
            rr.add_reactant(ecell4.core.Species(key), 1)
            rr.add_product(ecell4.core.Species(key), 1)
        else:
            raise RuntimeError(
                'unknown variable [{}] was used.'.format(key))
    exp = exp.format(*names)
    # print(exp)
    f = eval("lambda _r, _p, _v, _t, _rr: {0}".format(exp))
    f.__globals__.update(RATELAW_RESERVED_FUNCTIONS)
    f.__globals__.update((key, val) for key, val in RATELAW_RESERVED_CONSTANTS if val is not None)

    return (f, label)
    # return (lambda _r, _p, *args: eval(exp))

def parse_ReactionRule_options(elements):
    opts = {}

    if len(elements) == 0:
        raise RuntimeError(
            'only one attribute is allowed. [{:d}] given'.format(len(elements)))
    elif len(elements) == 1:
        opts['k'] = elements[0]
        return opts

    for elem in elements:
        if isinstance(elem, ecell4.core.ReactionRulePolicy):
            if 'policy' not in opts.keys():
                opts['policy'] = elem.get()
            else:
                opts['policy'] |= elem.get()
        else:
            if 'k' in opts.keys():
                raise RuntimeError('only one attribute is allowed. [%d] given' % (
                    len(elements)))
            opts['k'] = elem

    if 'k' not in opts.keys():
        raise RuntimeError('no kinetic rate or law is given.')

    return opts

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
        attribute_dimensionality = {'D': '[length]**2/[time]', 'radius': '[length]'}

        if not isinstance(obj, parseobj.OrExp):
            raise TypeError('An invalid object was given [{}]'.format(repr(obj)))

        elems = obj._elements()
        rhs = elems[-1]
        if isinstance(rhs, parseobj.ExpBase):
            return
        elif not isinstance(rhs, dict):
            raise TypeError('parameter must be given as a dict; "{}" given'.format(type(rhs)))

        for lhs in elems[: -1]:
            species_list = generate_Species(lhs)
            if len(species_list) != 1:
                raise ValueError(
                    'Only a single Species must be given; {:d} was given'.format(len(species_list)))
            elif species_list[0] is None:
                raise ValueError("No species given [{}]".format(repr(obj)))
            elif species_list[0][1] is not None:
                raise ValueError("Stoichiometry is not available here [{}]".format(repr(obj)))

            sp = species_list[0][0]

            for key, value in rhs.items():
                if unit.HAS_PINT and isinstance(value, unit._Quantity):
                    if (unit.STRICT and key in attribute_dimensionality
                        and not unit.check_dimensionality(value, attribute_dimensionality[key])):
                            raise ValueError("Cannot convert [{}] from '{}' ({}) to '{}'".format(
                                key, value.dimensionality, value.u, attribute_dimensionality[key]))
                    sp.set_attribute(key, value.to_base_units().magnitude)
                else:
                    sp.set_attribute(key, value)

            self.bitwise_operations.append(sp)

    def notify_comparisons(self, obj):
        raise RuntimeError(
            'ReactionRule definitions are not allowed'
            + ' in "species_attributes"')

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
        if not isinstance(obj, parseobj.CmpExp):
            raise RuntimeError('an invalid object was given [%s]' % (repr(obj)))
        elif isinstance(obj, parseobj.NeExp):
            warnings.warn('"<>" is deprecated; use "==" instead',
                          DeprecationWarning)

        lhs, rhs = obj._lhs, obj._rhs

        if isinstance(lhs, parseobj.OrExp):
            lhs = lhs._elements()[0]

        if not isinstance(rhs, parseobj.OrExp):
            raise RuntimeError('an invalid object was given'
                + ' as a right-hand-side [%s].' % (repr(rhs))
                + ' OrExp must be given')

        if len(rhs._elements()) == 0:
            raise RuntimeError('no product is given')

        opts = parse_ReactionRule_options(rhs._elements()[1: ])
        rhs = rhs._elements()[0]
        params = opts['k']

        lhs, rhs = generate_Species(lhs), generate_Species(rhs)
        lhs = tuple(sp for sp in lhs if sp is not None)
        rhs = tuple(sp for sp in rhs if sp is not None)

        if isinstance(obj, (parseobj.EqExp, parseobj.NeExp)):
            if not isinstance(params, (tuple, list)):
                raise RuntimeError(
                    'parameter must be a list or tuple with length 2; "%s" given'
                    % str(params))
            elif len(params) != 2:
                raise RuntimeError(
                    "parameter must be a list or tuple with length 2;"
                    + " length %d given" % len(params))
            # self.comparisons.append(generate_ReactionRule(lhs, rhs, params[0]))
            # self.comparisons.append(generate_ReactionRule(rhs, lhs, params[1]))
            rr = generate_ReactionRule(lhs, rhs, params[0])
            if 'policy' in opts.keys():
                rr.set_policy(opts['policy'])
            self.comparisons.append(rr)
            rr = generate_ReactionRule(rhs, lhs, params[1])
            if 'policy' in opts.keys():
                rr.set_policy(opts['policy'])
            self.comparisons.append(rr)
        elif isinstance(obj, parseobj.GtExp):
            # self.comparisons.append(generate_ReactionRule(lhs, rhs, params))
            rr = generate_ReactionRule(lhs, rhs, params)
            if 'policy' in opts.keys():
                rr.set_policy(opts['policy'])
            self.comparisons.append(rr)
        else:
            raise RuntimeError('an invalid object was given [%s]' % (repr(obj)))

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
    model : NetworkModel, NetfreeModel, or ODENetworkModel

    """
    try:
        if any([not isinstance(rr, ecell4.core.ReactionRule) for rr in REACTION_RULES]):
           from ecell4.ode import ODENetworkModel
           m = ODENetworkModel()
        elif seeds is not None or is_netfree:
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

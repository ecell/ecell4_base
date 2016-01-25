import copy
import types
import numbers
import warnings
import functools
import itertools

from . import parseobj
from .decorator_base import Callback, JustParseCallback, ParseDecorator

import ecell4.core

SEAMLESS_RATELAW_SUPPORT = False  #XXX: deprecated. use ENABLE_RATELAW
ENABLE_RATELAW = True
ENABLE_IMPLICIT_DECLARATION = True

PARAMETERS = []
SPECIES_ATTRIBUTES = []
REACTION_RULES = []

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

    if (callable(k)
        or ((SEAMLESS_RATELAW_SUPPORT or ENABLE_RATELAW) and isinstance(k, parseobj.ExpBase))
        or any([sp[1] is not None for sp in itertools.chain(lhs, rhs)])):
        from ecell4.ode import ODEReactionRule, ODERatelawCallback
        rr = ODEReactionRule()
        for sp in lhs:
            rr.add_reactant(sp[0], 1 if sp[1] is None else sp[1])
        for sp in rhs:
            rr.add_product(sp[0], 1 if sp[1] is None else sp[1])
        if (SEAMLESS_RATELAW_SUPPORT or ENABLE_RATELAW) and isinstance(k, parseobj.ExpBase):
            name = str(k)
            func = generate_ratelaw(k, rr)
            rr.set_ratelaw(ODERatelawCallback(func, name))
        elif callable(k):
            rr.set_ratelaw(ODERatelawCallback(k))
        else:
            rr.set_k(k)
        return rr
    elif isinstance(k, numbers.Number):
        return ecell4.core.ReactionRule([sp[0] for sp in lhs], [sp[0] for sp in rhs], k)

    raise RuntimeError(
        'parameter must be given as a number; "%s" given' % str(k))

def traverse_ParseObj(obj, keys):
    reserved_vars = ['_t', 'pi']
    reserved_funcs = ['exp', 'log', 'sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'pow']

    if isinstance(obj, parseobj.AnyCallable):
        obj = obj._as_ParseObj()

    if isinstance(obj, parseobj.ParseObj):
        if obj._size() == 1 and (
                obj._elems[0].name in reserved_funcs or
                obj._elems[0].name in reserved_vars):
            subobj = obj._elems[0]
            assert subobj.key is None
            assert subobj.modification is None
            if subobj.args is not None:
                assert subobj.name not in reserved_vars
                assert subobj.kwargs == {}
                subobj.args = tuple([
                    traverse_ParseObj(subobj.args[i], keys)
                    for i in range(len(subobj.args))])
            else:
                assert subobj.kwargs is None
        else:
            serial = ecell4.core.Species(str(obj)).serial()
            if serial in keys:
                return "{{{0:d}}}".format(keys.index(serial))
            keys.append(serial)
            return "{{{0:d}}}".format(len(keys) - 1)
    elif isinstance(obj, parseobj.ExpBase):
        for i in range(len(obj._elems)):
            obj._elems[i] = traverse_ParseObj(obj._elems[i], keys)
    return obj

def generate_ratelaw(obj, rr):
    keys = []
    exp = str(traverse_ParseObj(copy.deepcopy(obj), keys))
    aliases = {}
    for i, sp in enumerate(rr.reactants()):
        aliases[sp.serial()] = "_r[{0:d}]".format(i)
    for i, sp in enumerate(rr.products()):
        aliases[sp.serial()] = "_p[{0:d}]".format(i)
    names = []
    for key in keys:
        if key in aliases.keys():
            names.append(aliases[key])
        elif ENABLE_IMPLICIT_DECLARATION:
            names.append("_r[{0:d}]".format(len(rr.reactants())))
            aliases[key] = names[-1]
            rr.add_reactant(ecell4.core.Species(key), 1)
            rr.add_product(ecell4.core.Species(key), 1)
        else:
            raise RuntimeError(
                'unknown variable [{}] was used.'.format(key))
    exp = exp.format(*names)
    # print(exp)
    import math
    f = eval("lambda _r, _p, _v, _t, _rr: {0}".format(exp))
    f.__globals__['exp'] = math.exp
    f.__globals__['log'] = math.log
    f.__globals__['sin'] = math.sin
    f.__globals__['cos'] = math.cos
    f.__globals__['tan'] = math.tan
    f.__globals__['asin'] = math.asin
    f.__globals__['acos'] = math.acos
    f.__globals__['atan'] = math.atan
    f.__globals__['pi'] = math.pi
    f.__globals__['pow'] = pow
    return f
    # return (lambda _r, _p, *args: eval(exp))

class ParametersCallback(Callback):

    def __init__(self, *args):
        Callback.__init__(self)

        self.bitwise_operations = []

    def get(self):
        return copy.copy(self.bitwise_operations)

    def set(self):
        global PARAMETERS
        PARAMETERS.extend(self.bitwise_operations)

    def notify_bitwise_operations(self, obj):
        if not isinstance(obj, parseobj.OrExp):
            raise RuntimeError('an invalid object was given [%s]' % (repr(obj)))
        # elif len(obj._elements()) != 2:
        #     raise RuntimeError, 'only one attribute is allowed. [%d] given' % (
        #         len(obj._elements()))

        elems = obj._elements()
        rhs = elems[-1]
        if isinstance(rhs, parseobj.ExpBase):
            return

        for lhs in elems[: -1]:
            species_list = generate_Species(lhs)
            if len(species_list) != 1:
                raise RuntimeError(
                    'only a single species must be given; %d given'
                    % len(species_list))
            elif species_list[0] is None:
                raise RuntimeError("no species given [%s]" % (repr(obj)))
            elif species_list[0][1] is not None:
                raise RuntimeError(
                    "stoichiometry is not available here [%s]" % (repr(obj)))

            sp = species_list[0][0]

            if not isinstance(rhs, types.DictType):
                raise RuntimeError(
                    'parameter must be given as a dict; "%s" given'
                    % str(rhs))
            for key, value in rhs.items():
                if not (isinstance(key, types.StringType)
                    and isinstance(value, types.StringType)):
                    raise RuntimeError(
                        'attributes must be given as a pair of strings;'
                        + ' "%s" and "%s" given'
                        % (str(key), str(value)))
                sp.set_attribute(key, value)

            self.bitwise_operations.append(sp)

    def notify_comparisons(self, obj):
        raise RuntimeError(
            'ReactionRule definitions are not allowed'
            + ' in "world_inits"')

class SpeciesAttributesCallback(Callback):

    def __init__(self, *args):
        Callback.__init__(self)

        self.keys = None
        if len(args) > 0:
            for key in args:
                if not isinstance(key, (str, bytes)):
                    raise RuntimeError('non string key "%s" was given' % key)
            self.keys = args

        self.bitwise_operations = []

    def get(self):
        return copy.copy(self.bitwise_operations)

    def set(self):
        global SPECIES_ATTRIBUTES
        SPECIES_ATTRIBUTES.extend(self.bitwise_operations)

    def notify_bitwise_operations(self, obj):
        if not isinstance(obj, parseobj.OrExp):
            raise RuntimeError('an invalid object was given [%s]' % (repr(obj)))
        # elif len(obj._elements()) != 2:
        #     raise RuntimeError, 'only one attribute is allowed. [%d] given' % (
        #         len(obj._elements()))

        elems = obj._elements()
        rhs = elems[-1]
        if isinstance(rhs, parseobj.ExpBase):
            return

        for lhs in elems[: -1]:
            species_list = generate_Species(lhs)
            if len(species_list) != 1:
                raise RuntimeError(
                    'only a single species must be given; %d given'
                    % len(species_list))
            elif species_list[0] is None:
                raise RuntimeError("no species given [%s]" % (repr(obj)))
            elif species_list[0][1] is not None:
                raise RuntimeError(
                    "stoichiometry is not available here [%s]" % (repr(obj)))

            sp = species_list[0][0]

            if self.keys is None:
                if not isinstance(rhs, dict):
                    raise RuntimeError(
                        'parameter must be given as a dict; "%s" given'
                        % str(rhs))
                for key, value in rhs.items():
                    if not (isinstance(key, (str, bytes))
                        and isinstance(value, (str, bytes))):
                        raise RuntimeError(
                            'attributes must be given as a pair of strings;'
                            + ' "%s" and "%s" given'
                            % (str(key), str(value)))
                    sp.set_attribute(key, value)
            else:
                if not isinstance(rhs, (tuple, list)):
                    if len(self.keys) == 1:
                        rhs = (rhs, )
                    else:
                        raise RuntimeError(
                            'parameters must be given as a tuple or list; "%s" given'
                            % str(rhs))
                if len(rhs) != len(self.keys):
                    raise RuntimeError(
                        'the number of parameters must be %d; %d given'
                        % (len(self.keys), len(rhs)))
                else:
                    for key, value in zip(self.keys, rhs):
                        if not isinstance(value, (str, bytes)):
                            raise RuntimeError(
                                'paramter must be given as a string; "%s" given'
                                % str(value))
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
        elif len(rhs._elements()) != 2:
            raise RuntimeError('only one attribute is allowed. [%d] given' % (
                len(rhs._elements())))

        rhs, params = rhs._elements()
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
            self.comparisons.append(generate_ReactionRule(lhs, rhs, params[0]))
            self.comparisons.append(generate_ReactionRule(rhs, lhs, params[1]))
        elif isinstance(obj, parseobj.GtExp):
            self.comparisons.append(generate_ReactionRule(lhs, rhs, params))
        else:
            raise RuntimeError('an invalid object was given [%s]' % (repr(obj)))

def get_model(is_netfree=False, without_reset=False, seeds=None):
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

    Returns
    -------
    model : NetworkModel, NetfreeModel, or ODENetworkModel

    """
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
    for param in PARAMETERS:
        m.add_parameter(param)

    if not without_reset:
        reset_model()

    if seeds is not None:
        return m.expand(seeds)

    return m

def reset_model():
    """
    Reset all values, ``SPECIES_ATTRIBUTES`` and ``REACTIONRULES``,
    in the global scope.

    """
    global PARAMETERS
    global SPECIES_ATTRIBUTES
    global REACTION_RULES

    PARAMETERS = []
    SPECIES_ATTRIBUTES = []
    REACTION_RULES = []

reaction_rules = functools.partial(ParseDecorator, ReactionRulesCallback)
species_attributes = functools.partial(ParseDecorator, SpeciesAttributesCallback)
parameters = functools.partial(ParseDecorator, ParametersCallback)

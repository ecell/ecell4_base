import copy
import types
import numbers
import warnings
import functools
import itertools

from . import parseobj
from .decorator_base import Callback, JustParseCallback, parse_decorator, ParseDecorator

import ecell4.core


PARAMETERS = []
SPECIES_ATTRIBUTES = []
REACTION_RULES = []

def generate_Species(obj):
    if isinstance(obj, parseobj.AnyCallable):
        obj = obj._as_ParseObj()

    if isinstance(obj, parseobj.ParseObj):
        return (ecell4.core.Species(str(obj)), )
        # if len(elems) != 1:
        #     raise NotImplementedError, (
        #         'complex is not allowed yet; "%s"' % str(obj))
        # if (elems[0].args is not None
        #     or elems[0].kwargs is not None
        #     or elems[0].key is not None
        #     or elems[0].modification is not None):
        #     raise NotImplementedError, (
        #         'modification is not allowed yet; "%s"' % str(obj))
        # return (ecell4.core.Species(elems[0].name), )
    elif isinstance(obj, parseobj.InvExp):
        return (None, )
    elif isinstance(obj, parseobj.AddExp):
        subobjs = obj._elements()
        return tuple(itertools.chain(*[
            generate_Species(subobj) for subobj in subobjs]))
    else:
        raise RuntimeError('invalid expression; "%s" given' % str(obj))

def generate_ReactionRule(lhs, rhs, k=None):
    if k is None:
        raise RuntimeError('no parameter is specified')
    elif callable(k):
        from ecell4.ode import ODEReactionRule, ODERatelawCallback
        rr = ODEReactionRule()
        for sp in lhs:
            rr.add_reactant(sp, 1)
        for sp in rhs:
            rr.add_product(sp, 1)
        rr.set_ratelaw(ODERatelawCallback(k))
        return rr
    elif isinstance(params, numbers.Number):
        return ecell4.core.ReactionRule(lhs, rhs, k)

    raise RuntimeError(
        'parameter must be given as a number; "%s" given'
        % str(params))

# def generate_ReactionRule(lhs, rhs, k=None):
#     return ecell4.core.ReactionRule(lhs, rhs, k)
#     # if len(lhs) == 0:
#     #     if len(rhs) != 1:
#     #         raise RuntimeError(
#     #             "the number of products must be 1; %d given" % len(rhs))
#     #     return ecell4.core.create_synthesis_reaction_rule(rhs[0], k)
#     # elif len(lhs) == 1:
#     #     if len(rhs) == 0:
#     #         return ecell4.core.create_degradation_reaction_rule(lhs[0], k)
#     #     elif len(rhs) == 1:
#     #         return ecell4.core.create_unimolecular_reaction_rule(
#     #             lhs[0], rhs[0], k)
#     #     elif len(rhs) == 2:
#     #         return ecell4.core.create_unbinding_reaction_rule(
#     #             lhs[0], rhs[0], rhs[1], k)
#     #     else:
#     #         raise RuntimeError(
#     #             "the number of products must be less than 3; %d given"
#     #             % len(rhs))
#     # elif len(lhs) == 2:
#     #     if len(rhs) == 1:
#     #         return ecell4.core.create_binding_reaction_rule(
#     #             lhs[0], lhs[1], rhs[0], k)
#     #     else:
#     #         raise RuntimeError(
#     #             "the number of products must be 1; %d given" % len(rhs))
#     # raise RuntimeError(
#     #     "the number of reactants must be less than 3; %d given" % len(lhs))

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

            sp = species_list[0]
            if sp is None:
                raise RuntimeError("no species given [%s]" % (repr(obj)))

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

            sp = species_list[0]
            if sp is None:
                raise RuntimeError("no species given [%s]" % (repr(obj)))

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
    global PARAMETERS
    global SPECIES_ATTRIBUTES
    global REACTION_RULES

    PARAMETERS = []
    SPECIES_ATTRIBUTES = []
    REACTION_RULES = []

reaction_rules = functools.partial(ParseDecorator, ReactionRulesCallback)
species_attributes = functools.partial(ParseDecorator, SpeciesAttributesCallback)
parameters = functools.partial(ParseDecorator, ParametersCallback)

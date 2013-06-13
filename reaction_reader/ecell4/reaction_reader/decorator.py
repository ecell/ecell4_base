import copy
import types
import numbers
import warnings
import functools

import parseobj

import ecell4.core


def generate_Species(obj):
    if isinstance(obj, parseobj.AnyCallable):
        obj = obj._as_ParseObj()

    if isinstance(obj, parseobj.ParseObj):
        elems = obj._get_elements()
        if len(elems) != 1:
            raise NotImplementedError, (
                'complex is not allowed yet; "%s"' % str(obj))
        if (elems[0].args is not None
            or elems[0].kwargs is not None
            or elems[0].key is not None
            or elems[0].modification is not None):
            raise NotImplementedError, (
                'modification is not allowed yet; "%s"' % str(obj))
        if elems[0].inv:
            return ((None, elems[0].param), )
        else:
            return ((ecell4.core.Species(elems[0].name), elems[0].param), )
    elif isinstance(obj, parseobj.ParseObjSet):
        subobjs = obj._get_objects()
        return tuple(generate_Species(subobj)[0] for subobj in subobjs)
    raise RuntimeError, 'invalid expression; "%s" given' % str(obj)

def generate_ReactionRule(lhs, rhs, k=0.0):
    if len(lhs) == 0:
        if len(rhs) != 1:
            raise RuntimeError, (
                "the number of products must be 1; %d given" % len(rhs))
        return ecell4.core.create_synthesis_reaction_rule(rhs[0], k)
    elif len(lhs) == 1:
        if len(rhs) == 0:
            return ecell4.core.create_degradation_reaction_rule(lhs[0], k)
        elif len(rhs) == 1:
            return ecell4.core.create_unimolecular_reaction_rule(
                lhs[0], rhs[0], k)
        elif len(rhs) == 2:
            return ecell4.core.create_unbinding_reaction_rule(
                lhs[0], rhs[0], rhs[1], k)
        else:
            raise RuntimeError, (
                "the number of products must be less than 3; %d given"
                % len(rhs))
    elif len(lhs) == 2:
        if len(rhs) == 1:
            return ecell4.core.create_binding_reaction_rule(
                lhs[0], lhs[1], rhs[0], k)
        else:
            raise RuntimeError, (
                "the number of products must be 1; %d given" % len(rhs))
    raise RuntimeError, (
        "the number of reactants must be less than 3; %d given" % len(lhs))

class SpeciesAttributesCallback(object):

    def __init__(self, *args):
        self.keys = None
        if len(args) > 0:
            for key in args:
                if not isinstance(key, types.StringType):
                    raise RuntimeError, 'non string key "%s" was given' % key
            self.keys = args

        self.bitwise_operations = []

    def get(self):
        return copy.copy(self.bitwise_operations)

    def notify_unary_operations(self, optr, target):
        pass

    def notify_bitwise_operations(self, optr, lhs, rhs):
        if optr != '|':
            raise RuntimeError, 'operator "%s" not allowed' % optr

        species_list = generate_Species(lhs)
        if len(species_list) != 1:
            raise RuntimeError, (
                'only a single species must be given; %d given'
                % len(species_list))

        sp, _ = species_list[0]
        if sp is None:
            raise RuntimeError, 'never use "~" in "species_attributes"'

        if self.keys is None:
            if not isinstance(rhs, types.DictType):
                raise RuntimeError, (
                    'parameter must be given as a dict; "%s" given'
                    % str(rhs))
            for key, value in rhs.items():
                if not (isinstance(key, types.StringType)
                    and isinstance(value, types.StringType)):
                    raise RuntimeError, (
                        'attributes must be given as a pair of strings;'
                        + ' "%s" and "%s" given'
                        % (str(key), str(value)))
                sp.set_attribute(key, value)
        else:
            if not (isinstance(rhs, types.TupleType)
                and isinstance(rhs, types.ListType)):
                if len(self.keys) == 1:
                    rhs = (rhs, )
                else:
                    raise RuntimeError, (
                        'parameters must be given as a tuple or list; "%s" given'
                        % str(rhs))
            if len(rhs) != len(self.keys):
                raise RuntimeError, (
                    'the number of parameters must be %d; %d given'
                    % (len(self.keys), len(rhs)))
            else:
                for key, value in zip(self.keys, rhs):
                    if not isinstance(value, types.StringType):
                        raise RuntimeError, (
                            'paramter must be given as a string; "%s" given'
                            % str(value))
                    sp.set_attribute(key, value)

        self.bitwise_operations.append(sp)

    def notify_comparisons(self, optr, lhs, rhs):
        raise RuntimeError, (
            'ReactionRule definitions are not allowed'
            + ' in "species_attributes"')

class ReactionRulesCallback(object):

    def __init__(self):
        self.comparisons = []

    def get(self):
        return copy.copy(self.comparisons)

    def notify_unary_operations(self, optr, target):
        pass

    def notify_bitwise_operations(self, optr, lhs, rhs):
        pass

    def notify_comparisons(self, optr, lhs, rhs):
        if optr == "!=":
            warnings.warn('"<>" is deprecated; use "==" instead',
                          DeprecationWarning)

        lhs, rhs = generate_Species(lhs), generate_Species(rhs)
        lhs = tuple(sp for sp, param in lhs if sp is not None)
        k = rhs[-1][1]
        rhs = tuple(sp for sp, param in rhs if sp is not None)

        if optr == "!=" or optr == "==":
            if not (isinstance(k, types.ListType)
                or isinstance(k, types.TupleType)):
                raise RuntimeError, (
                    'parameter must be a list or tuple with length 2; "%s" given'
                    % str(k))
            elif len(k) != 2:
                raise RuntimeError, (
                    "parameter must be a list or tuple with length 2;"
                    + " length %d given" % len(k))
            elif not (isinstance(k[0], numbers.Number)
                and isinstance(k[1], numbers.Number)):
                raise RuntimeError, (
                    'parameters must be given as a list or tuple of numbers;'
                    + ' "%s" given' % str(k))
            self.comparisons.append(generate_ReactionRule(lhs, rhs, k[0]))
            self.comparisons.append(generate_ReactionRule(rhs, lhs, k[1]))
        elif optr == ">":
            if k is None:
                raise RuntimeError, 'no parameter is specified'
            elif not isinstance(k, numbers.Number):
                raise RuntimeError, (
                    'parameter must be given as a number; "%s" given' % str(k))
            self.comparisons.append(generate_ReactionRule(lhs, rhs, k))
        else:
            raise RuntimeError, 'operator "%s" not allowed' % optr

class JustParseCallback(object):

    def __init__(self):
        self.comparisons = []
        # self.bitwise_optrs = []

    def get(self):
        return copy.copy(self.comparisons)

    def notify_unary_operations(self, optr, target):
        pass

    def notify_bitwise_operations(self, optr, lhs, rhs):
        # self.bitwise_optrs.append((optr, lhs, rhs))
        pass

    def notify_comparisons(self, optr, lhs, rhs):
        self.comparisons.append((optr, lhs, rhs))

class Callback(object):
    """callback before the operations"""

    def __init__(self):
        # self.bitwise_operations = []
        self.comparisons = []

    def get(self):
        return copy.copy(self.comparisons)

    def notify_unary_operations(self, optr, target):
        pass

    def notify_bitwise_operations(self, optr, lhs, rhs):
        # self.bitwise_operations.append((optr, lhs, rhs))
        pass

    def notify_comparisons(self, optr, lhs, rhs):
        if optr == "!=":
            warnings.warn('"<>" is deprecated; use "==" instead',
                          DeprecationWarning)
        self.comparisons.append((optr, lhs, rhs))

def parse_decorator(callback_class, func):
    @functools.wraps(func)
    def wrapped(*args, **kwargs):
        cache = callback_class()
        vardict = copy.copy(func.func_globals)
        for k in func.func_code.co_names:
            if (not k in vardict.keys()
                and not k in dir(vardict['__builtins__'])): # is this enough?
                vardict[k] = parseobj.AnyCallable(cache, k)
        g = types.FunctionType(func.func_code, vardict)
        with warnings.catch_warnings():
            # warnings.simplefilter("always")
            g(*args, **kwargs)
        return cache.get()
    return wrapped

# reaction_rules = functools.partial(parse_decorator, Callback)
reaction_rules = functools.partial(parse_decorator, ReactionRulesCallback)
species_attributes = functools.partial(parse_decorator, SpeciesAttributesCallback)
just_parse = functools.partial(parse_decorator, JustParseCallback)

def species_attributes_with_keys(*args):
    def create_callback():
        return SpeciesAttributesCallback(*args)
    return functools.partial(parse_decorator, create_callback)

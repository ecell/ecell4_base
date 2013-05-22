import copy
import types
import numbers
import warnings
from functools import wraps

import parseobj

import ecell4.core


def generate_Species(obj):
    if isinstance(obj, parseobj.AnyCallable):
        obj = obj._as_ParseObj()

    if isinstance(obj, parseobj.ParseObj):
        elems = obj._get_elements()
        if len(elems) != 1:
            raise RuntimeError
        if (elems[0].args is not None
            or elems[0].kwargs is not None
            or elems[0].key is not None):
            raise RuntimeError
        return ((ecell4.core.Species(elems[0].name), elems[0].param), )
    elif isinstance(obj, parseobj.ParseObjSet):
        subobjs = obj._get_objects()
        return tuple(generate_Species(subobj)[0] for subobj in subobjs)
    raise RuntimeError

def generate_ReactionRule(lhs, rhs, k=0.0):
    if len(lhs) == 0:
        if len(rhs) != 1:
            raise RuntimeError
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
            raise RuntimeError
    elif len(lhs) == 2:
        if len(rhs) == 1:
            return ecell4.core.create_binding_reaction_rule(
                lhs[0], lhs[1], rhs[0], k)
        else:
            raise RuntimeError
    raise RuntimeError

class ReactionRuleCallback(object):

    def __init__(self):
        self.comparisons = []

    def notify_bitwise_operations(self, optr, lhs, rhs):
        pass

    def notify_comparisons(self, optr, lhs, rhs):
        if optr == "!=":
            warnings.warn('"<>" is deprecated; use "==" instead',
                          DeprecationWarning)

        lhs, rhs = generate_Species(lhs), generate_Species(rhs)
        lhs = tuple(sp for sp, param in lhs)
        k = rhs[-1][1]
        rhs = tuple(sp for sp, param in rhs)

        if optr == "!=" or optr == "==":
            if not (isinstance(k, types.ListType)
                or isinstance(k, types.TupleType)):
                raise RuntimeError
            elif len(k) != 2:
                raise RuntimeError
            elif not (isinstance(k[0], numbers.Number)
                and isinstance(k[1], numbers.Number)):
                raise RuntimeError
            self.comparisons.append(generate_ReactionRule(lhs, rhs, k[0]))
            self.comparisons.append(generate_ReactionRule(rhs, lhs, k[1]))
        else:
            if not isinstance(k, numbers.Number):
                raise RuntimeError
            self.comparisons.append(generate_ReactionRule(lhs, rhs, k))

class Callback(object):
    """callback before the operations"""

    def __init__(self):
        # self.bitwise_operations = []
        self.comparisons = []

    def notify_bitwise_operations(self, optr, lhs, rhs):
        # self.bitwise_operations.append((optr, lhs, rhs))
        pass

    def notify_comparisons(self, optr, lhs, rhs):
        if optr == "!=":
            warnings.warn('"<>" is deprecated; use "==" instead',
                          DeprecationWarning)
        self.comparisons.append((optr, lhs, rhs))

def reaction_rules(func):
    @wraps(func)
    def wrapped(*args, **kwargs):
        cache = Callback()
        # cache = ReactionRuleCallback()
        vardict = copy.copy(globals())
        for k in func.func_code.co_names:
            if not k in vardict.keys():
                vardict[k] = parseobj.AnyCallable(cache, k)
        g = types.FunctionType(func.func_code, vardict)
        with warnings.catch_warnings():
            # warnings.simplefilter("always")
            g(*args, **kwargs)
        return cache.comparisons
    return wrapped

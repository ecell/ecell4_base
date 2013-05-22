import copy
import types
import warnings
from functools import wraps

from parseobj import AnyCallable


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
        vardict = copy.copy(globals())
        for k in func.func_code.co_names:
            if not k in vardict.keys():
                vardict[k] = AnyCallable(cache, k)
        g = types.FunctionType(func.func_code, vardict)
        with warnings.catch_warnings():
            # warnings.simplefilter("always")
            g(*args, **kwargs)
        return cache.comparisons
    return wrapped

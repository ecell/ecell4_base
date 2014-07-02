import copy
import types
import warnings
import functools

import parseobj


class Callback(object):
    """callback before the operations"""

    def __init__(self):
        pass

    def get(self):
        return None

    def notify_unary_operations(self, obj):
        pass

    def notify_bitwise_operations(self, obj):
        pass

    def notify_comparisons(self, obj):
        pass

class JustParseCallback(Callback):

    def __init__(self):
        Callback.__init__(self)

        self.comparisons = []

    def get(self):
        return copy.copy(self.comparisons)

    def notify_comparisons(self, obj):
        if isinstance(obj, parseobj.NeExp):
            warnings.warn('"<>" is deprecated; use "==" instead',
                          DeprecationWarning)
        self.comparisons.append(obj)

def parse_decorator(callback_class, func):
    @functools.wraps(func)
    def wrapped(*args, **kwargs):
        cache = callback_class()
        vardict = copy.copy(func.func_globals)
        for ignore in ("_", "__", "___", "_i", "_ii", "_iii",
            "_i1", "_i2", "_i3", "_dh", "_sh", "_oh"):
            if ignore in vardict.keys():
                del vardict[ignore]
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

just_parse = functools.partial(parse_decorator, JustParseCallback)

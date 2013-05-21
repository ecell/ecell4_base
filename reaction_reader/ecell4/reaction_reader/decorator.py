import copy
import types
from functools import wraps

from parseobj import AnyCallable


def reaction_rules(func):
    @wraps(func)
    def wrapped(*args, **kwargs):
        cache = []
        vardict = copy.copy(globals())
        for k in func.func_code.co_names:
            if not k in vardict.keys():
                vardict[k] = AnyCallable(cache, k)
        g = types.FunctionType(func.func_code, vardict)
        # return [rr for rr in g(*args, **kwargs)]
        g(*args, **kwargs)
        return cache
    return wrapped

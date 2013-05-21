import copy
import types
from parseobj import AnyCallable


def reaction_rules(func):
    def wrapper(*args, **kwargs):
        cache = []
        vardict = copy.copy(globals())
        for k in func.func_code.co_names:
            if not k in vardict.keys():
                vardict[k] = AnyCallable(cache, k)
        g = types.FunctionType(func.func_code, vardict)
        # return [rr for rr in g(*args, **kwargs)]
        g(*args, **kwargs)
        return cache
    return wrapper

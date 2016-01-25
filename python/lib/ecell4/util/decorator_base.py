import copy
import types
import warnings
import functools
import inspect

from . import parseobj


class Callback(object):
    """callback before the operations"""

    def __init__(self):
        pass

    def set(self):
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

def keys_from_builtins(vardict):
    b = vardict['__builtins__']
    if isinstance(b, types.ModuleType):
        return dir(b)
    else:
        return dict(b).keys()

class ParseDecorator:

    def __init__(self, callback_class, func=None):
        self.__callback_class = callback_class

        self.__callback = None
        self.__newvars = {}

        if func is not None:
            self.__func = func
            functools.update_wrapper(self, func)
        else:
            self.__func = lambda *args, **kwargs: []

    def __call__(self, *args, **kwargs):
        self.__callback = self.__callback_class()
        try:
            vardict = copy.copy(self.__func.func_globals)
            func_code = self.__func.func_code
            name = self.__func.func_name
            defaults = self.__func.func_defaults
        except AttributeError:
            vardict = copy.copy(self.__func.__globals__)
            func_code = self.__func.__code__
            name = self.__func.__name__
            defaults = self.__func.__defaults__

        ignores = ("_", "__", "___", "_i", "_ii", "_iii",
            "_i1", "_i2", "_i3", "_dh", "_sh", "_oh")
        for ignore in ignores:
            if ignore in vardict.keys():
                del vardict[ignore]
        for k in func_code.co_names:
            if k == '_eval':
                vardict[k] = self.__evaluate
            elif (not k in vardict.keys()
                and not k in keys_from_builtins(vardict)): # is this enough?
                vardict[k] = parseobj.AnyCallable(self.__callback, k)
        g = types.FunctionType(func_code, vardict, name=name, argdefs=defaults)
        with warnings.catch_warnings():
            # warnings.simplefilter("always")
            g(*args, **kwargs)
        return self.__callback.get()

    def __enter__(self):
        # print "ParseDecorator#__enter__"
        self.__callback = self.__callback_class()
        calling_frame = inspect.currentframe().f_back
        vardict = copy.copy(calling_frame.f_globals)
        ignores = ("_", "__", "___", "_i", "_ii", "_iii",
            "_i1", "_i2", "_i3", "_dh", "_sh", "_oh")
        for k in calling_frame.f_code.co_names:
            if k == '_eval':
                calling_frame.f_globals[k] = self.__evaluate
                self.__newvars[k] = None
            elif k in ignores:
                # print "WARNING: '%s' was overridden." % k
                calling_frame.f_globals[k] = parseobj.AnyCallable(self.__callback, k)
                self.__newvars[k] = vardict.get(k)
            elif (not k in vardict.keys()
                and not k in keys_from_builtins(vardict)):
                # print "WARNING: '%s' is undefined." % k
                calling_frame.f_globals[k] = parseobj.AnyCallable(self.__callback, k)
                self.__newvars[k] = None
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        # print "ParseDecorator#__exit__", exc_type, exc_value, traceback
        if self.__callback is not None:
            if exc_type is None:
                # print self.__callback.get()
                self.__callback.set()

            self.__callback = None
            calling_frame = inspect.currentframe().f_back
            for k, v in self.__newvars.items():
                if v is None:
                    del calling_frame.f_globals[k]
                    # print "WARNING: '%s' was removed." % k
                else:
                    calling_frame.f_globals[k] = v
                    # print "WARNING: '%s' was recovered to be '%s'." % (k, v)

    def __evaluate(self, expr, params={}):
        class AnyCallableLocals:

            def __init__(self, callback, locals):
                self.callback = callback
                self.locals = locals

            def __getitem__(self, key):
                if key in self.locals.keys():
                    return self.locals[key]
                return parseobj.AnyCallable(self.callback, key)
        l = locals()
        l.update(params)
        return eval(expr, globals(), AnyCallableLocals(self.__callback, l))

def parse_decorator(callback_class, func):
    @functools.wraps(func)
    def wrapped(*args, **kwargs):
        cache = callback_class()
        try:
            vardict = copy.copy(self.__func.func_globals)
            func_code = func.func_code
            name = self.__func.func_name
            defaults = self.__func.func_defaults
        except AttributeError:
            vardict = copy.copy(self.__func.__globals__)
            func_code = func.__code__
            name = self.__func.__name__
            defaults = self.__func.__defaults__
        for ignore in ("_", "__", "___", "_i", "_ii", "_iii",
            "_i1", "_i2", "_i3", "_dh", "_sh", "_oh"):
            if ignore in vardict.keys():
                del vardict[ignore]
        for k in func_code.co_names:
            if (not k in vardict.keys()
                and not k in keys_from_builtins(vardict)): # is this enough?
                vardict[k] = parseobj.AnyCallable(cache, k)
        g = types.FunctionType(func_code, vardict, name=name, argdefs=defaults)
        with warnings.catch_warnings():
            # warnings.simplefilter("always")
            g(*args, **kwargs)
        return cache.get()
    return wrapped

# just_parse = functools.partial(parse_decorator, JustParseCallback)
just_parse = functools.partial(ParseDecorator, JustParseCallback)

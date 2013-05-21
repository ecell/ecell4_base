import operator
from logger import log_call


class AnyCallable:
    """AnyCallable must be immutable."""

    def __init__(self, root, name):
        self.__root = root # a reference to cache
        self.__name = name

    def as_ParseObj(self):
        return ParseObj(self.__root, self.__name)

    def __call__(self, *args, **kwargs):
        retval = self.as_ParseObj()
        return retval(*args, **kwargs)

    def __getitem__(self, key):
        retval = self.as_ParseObj()
        return retval[key]

    def __getattr__(self, key):
        return getattr(self.as_ParseObj(), key)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return str(self.as_ParseObj())

class ParseElem:

    def __init__(self, name):
        self.name = name
        self.args = []
        self.kwargs = []
        self.key = None
        self.param = None

    @log_call
    def __call__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        return self

    def __getitem__(self, key):
        self.key = key
        return self

    def __coerce__(self, other):
        return None

    @log_call
    def __or__(self, rhs):
        self.param = rhs
        return self

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        label = self.name
        if len(self.args) > 0 or len(self.kwargs) > 0:
            attrs = ["%s" % v for v in sorted(self.args)]
            attrs += ["%s=%s" % (k, self.kwargs[k]) for k in sorted(self.kwargs.keys())]
            label += "(%s)" % (",".join(attrs))
        if self.key is not None:
            label += "[%s]" % self.key
        if self.param is not None:
            label += "|%s" % self.param
        return label

class ParseObj:
    """All the members in ParseObj must start from "__"."""

    def __init__(self, root, name, elems=[]):
        self.__root = root # a reference to cache
        self.__elems = elems + [ParseElem(name)]

    @log_call
    def __call__(self, *args, **kwargs):
        self.__elems[-1] = self.__elems[-1](*args, **kwargs)
        return self

    @log_call
    def __add__(self, rhs):
        return ParseObjSet(self.__root, (self, rhs))

    def __getitem__(self, key):
        self.__elems[-1] = self.__elems[-1][key]
        return self

    @log_call
    def __getattr__(self, key):
        return ParseObjPartial(self, key)

    def __coerce__(self, other):
        return None

    @log_call
    def __or__(self, rhs):
        self.__elems[-1] = (self.__elems[-1] | rhs)
        return self

    @log_call
    def __gt__(self, rhs):
        self.__root.append((self, rhs))
        return (self, rhs)

    @log_call
    def __lshift__(self, other):
        """not for users, but only for ParseObjPartial"""
        self.__elems.append(other)
        return self

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        # XXX: sort by element's name
        # XXX: separate params
        return ".".join([str(elem) for elem in self.__elems])

class ParseObjSet:

    def __init__(self, root, objs):
        if len(objs) < 2:
            raise RuntimeError
        self.__root = root
        self.__objs = list(objs)

    def __coerce__(self, other):
        return None

    @log_call
    def __or__(self, rhs):
        self.__objs[-1] = (self.__objs[-1] | rhs)
        return self

    @log_call
    def __gt__(self, rhs):
        self.__root.append((self, rhs))
        return (self, rhs)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        labels = [str(obj) for obj in self.__objs]
        return "+".join(labels)

class ParseObjPartial:

    def __init__(self, lhs, name):
        self.__lhs = lhs
        self.__name = name

    @log_call
    def __call__(self, *args, **kwargs):
        rhs = ParseElem(self.__name)
        rhs = rhs(*args, **kwargs)
        return operator.lshift(self.__lhs, rhs) # (self.lhs << rhs)

from logger import log_call
import operator


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
        self.args = None
        self.kwargs = None
        self.key = None
        self.param = None

    @log_call
    def set_arguments(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def set_key(self, key):
        self.key = key

    @log_call
    def set_parameter(self, rhs):
        self.param = rhs

    def __coerce__(self, other):
        return None

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        label = self.name

        if self.args is not None or self.kwargs is not None:
            attrs = []
            if self.args is not None:
                attrs += ["%s" % str(v) for v in self.args]
            if self.kwargs is not None:
                attrs += ["%s=%s" % (k, v) for k, v in self.kwargs.items()]
            label += "(%s)" % (",".join(attrs))

        if self.key is not None:
            label += "[%s]" % str(self.key)
        if self.param is not None:
            label += "|%s" % str(self.param)
        return label

class ParseObj:
    """All the members in ParseObj must start from "__"."""

    def __init__(self, root, name, elems=[]):
        self.__root = root # a reference to cache
        self.__elems = elems + [ParseElem(name)]

    @log_call
    def __call__(self, *args, **kwargs):
        self.__elems[-1].set_arguments(*args, **kwargs)
        return self

    @log_call
    def __add__(self, rhs):
        if isinstance(rhs, ParseObj):
            return ParseObjSet(self.__root, (self, rhs))
        elif isinstance(rhs, ParseObjSet):
            return ParseObjSet(self.__root, operator.rshift(rhs, [self]))
        raise RuntimeError, "never get here"

    def __getitem__(self, key):
        self.__elems[-1].set_key(key)
        return self

    @log_call
    def __getattr__(self, key):
        return ParseObjPartial(self, key)

    @log_call
    def __or__(self, rhs):
        self.__elems[-1].set_parameter(rhs)
        return self

    @log_call
    def __gt__(self, rhs):
        retval = ("gt", self, rhs)
        self.__root.append(retval)
        return retval

    @log_call
    def __eq__(self, rhs):
        retval = ("eq", self, rhs)
        self.__root.append(retval)
        return retval

    @log_call
    def __ne__(self, rhs):
        retval = ("neq", self, rhs)
        self.__root.append(retval)
        return retval

    @log_call
    def __lshift__(self, other):
        """not for users, but only for ParseObjPartial"""
        self.__elems.append(other)
        return self

    def __rshift__(self, other):
        if type(other) is not list:
            raise RuntimeError
        other.extend(self.__elems)
        return other

    def __coerce__(self, other):
        return None

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        labels = [str(elem) for elem in self.__elems]
        return ".".join(labels)

class ParseObjSet:

    def __init__(self, root, objs):
        if len(objs) < 2:
            raise RuntimeError
        self.__root = root
        self.__objs = list(objs)

    def __call__(self, *args, **kwargs):
        raise RuntimeError

    def __getitem__(self, key):
        raise RuntimeError

    def __getattr__(self, key):
        raise RuntimeError

    def __lshift__(self, other):
        raise RuntimeError

    @log_call
    def __add__(self, rhs):
        if isinstance(rhs, ParseObj):
            self.__objs.append(rhs)
            return self
        elif isinstance(rhs, ParseObjSet):
            operator.rshift(rhs, self.__objs)
            return self
        raise RuntimeError, "never get here"

    @log_call
    def __or__(self, rhs):
        self.__objs[-1] = (self.__objs[-1] | rhs)
        return self

    @log_call
    def __gt__(self, rhs):
        retval = ("gt", self, rhs)
        self.__root.append(retval)
        return retval

    @log_call
    def __eq__(self, rhs):
        retval = ("eq", self, rhs)
        self.__root.append(retval)
        return retval

    @log_call
    def __ne__(self, rhs):
        retval = ("neq", self, rhs)
        self.__root.append(retval)
        return retval

    def __rshift__(self, other):
        if type(other) is not list:
            raise RuntimeError
        for obj in self.__objs:
            other.extend(operator.rshift(obj, list()))
        return other

    def __coerce__(self, other):
        return None

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
        rhs.set_arguments(*args, **kwargs)
        return operator.lshift(self.__lhs, rhs) # (self.lhs << rhs)

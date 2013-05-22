from logger import log_call
import operator


class AnyCallable:
    """AnyCallable must be immutable.
All the members must start with '_'."""

    def __init__(self, root, name):
        self.__root = root # a reference to cache
        self.__name = name

    def _as_ParseObj(self):
        return ParseObj(self.__root, self.__name)

    def __getattr__(self, key):
        return getattr(self._as_ParseObj(), key)

    def __coerce__(self, other):
        return None

    def __str__(self):
        return self.__name

    def __repr__(self):
        return "<%s.%s: %s>" % (
            self.__class__.__module__, self.__class__.__name__, str(self))

class ParseElem:

    def __init__(self, name):
        self.name = name
        self.args = None
        self.kwargs = None
        self.key = None
        self.param = None

    def set_arguments(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def set_key(self, key):
        self.key = key

    def set_parameter(self, rhs):
        self.param = rhs

    def __str__(self):
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

    def __repr__(self):
        return "<%s.%s: %s>" % (
            self.__class__.__module__, self.__class__.__name__, str(self))

class ParseObj:
    """All the members must start with '_'."""

    def __init__(self, root, name, elems=[]):
        self.__root = root # a reference to cache
        self.__elems = elems + [ParseElem(name)]

    @log_call
    def __call__(self, *args, **kwargs):
        self.__elems[-1].set_arguments(*args, **kwargs)
        return self

    @log_call
    def __add__(self, rhs):
        if isinstance(rhs, AnyCallable):
            return ParseObjSet(self.__root, (self, rhs._as_ParseObj()))
        elif isinstance(rhs, ParseObj):
            return ParseObjSet(self.__root, (self, rhs))
        elif isinstance(rhs, ParseObjSet):
            return ParseObjSet(self.__root, rhs._get_objects() + [self])
        raise RuntimeError, "never get here"

    @log_call
    def __getitem__(self, key):
        self.__elems[-1].set_key(key)
        return self

    @log_call
    def __getattr__(self, key):
        if key[0] == "_":
            # raise AttributeError, (
            #     "'%s' object has no attribute '%s'"
            #         % (self.__class__.__name__, key))
            raise RuntimeError, (
                "'%s' object has no attribute '%s'"
                    % (self.__class__.__name__, key))
        self.__elems.append(ParseElem(key))
        return self

    @log_call
    def __or__(self, rhs):
        optr = "|"
        self.__root.notify_bitwise_operations(optr, self, rhs)

        self.__elems[-1].set_parameter(rhs)
        return self

    @log_call
    def __gt__(self, rhs):
        optr = ">"
        self.__root.notify_comparisons(optr, self, rhs)
        return (optr, self, rhs)

    @log_call
    def __eq__(self, rhs):
        optr = "=="
        self.__root.notify_comparisons(optr, self, rhs)
        return (optr, self, rhs)

    @log_call
    def __ne__(self, rhs):
        optr = "!="
        self.__root.notify_comparisons(optr, self, rhs)
        return (optr, self, rhs)

    def __coerce__(self, other):
        return None

    def __str__(self):
        labels = [str(elem) for elem in self.__elems]
        return ".".join(labels)

    def __repr__(self):
        return "<%s.%s: %s>" % (
            self.__class__.__module__, self.__class__.__name__, str(self))

class ParseObjSet:
    """All the members must start with '_'."""

    def __init__(self, root, objs):
        if len(objs) < 2:
            raise RuntimeError
        self.__root = root
        self.__objs = list(objs)

    def _get_objects(self):
        return copy.copy(self.__objs)

    def __call__(self, *args, **kwargs):
        raise RuntimeError

    def __getitem__(self, key):
        raise RuntimeError

    def __getattr__(self, key):
        raise RuntimeError

    @log_call
    def __add__(self, rhs):
        if isinstance(rhs, AnyCallable):
            self.__objs.append(rhs._as_ParseObj())
            return self
        elif isinstance(rhs, ParseObj):
            self.__objs.append(rhs)
            return self
        elif isinstance(rhs, ParseObjSet):
            self.__objs.extend(rhs._get_objects())
            return self
        raise RuntimeError, "never get here"

    @log_call
    def __or__(self, rhs):
        optr = "|"
        self.__root.notify_bitwise_operations(optr, self, rhs)

        self.__objs[-1] = (self.__objs[-1] | rhs)
        return self

    @log_call
    def __gt__(self, rhs):
        optr = ">"
        self.__root.notify_comparisons(optr, self, rhs)
        return (optr, self, rhs)

    @log_call
    def __eq__(self, rhs):
        optr = "=="
        self.__root.notify_comparisons(optr, self, rhs)
        return (optr, self, rhs)

    @log_call
    def __ne__(self, rhs):
        optr = "!="
        self.__root.notify_comparisons(optr, self, rhs)
        return (optr, self, rhs)

    def __coerce__(self, other):
        return None

    def __str__(self):
        labels = [str(obj) for obj in self.__objs]
        return "+".join(labels)

    def __repr__(self):
        return "<%s.%s: %s>" % (
            self.__class__.__module__, self.__class__.__name__, str(self))

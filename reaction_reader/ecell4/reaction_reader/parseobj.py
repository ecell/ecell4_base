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

class ParseObj:

    def __init__(self, root, name, lhs=None):
        self.root = root # a reference to cache

        self.lhs = lhs
        self.name = name
        self.args = []
        self.kwargs = []
        self.key = None
        self.param = None
        self.called = False

    @log_call
    def __call__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        self.called = True
        return self

    @log_call
    def __add__(self, rhs):
        return ParseObjSet((self, rhs))

    def __getitem__(self, key):
        if self.called:
            raise RuntimeError
        self.key = key
        return self

    def __getattr__(self, key):
        return ParseObjPartial(self, key)

    def __coerce__(self, other):
        return None

    @log_call
    def __or__(self, rhs):
        self.param = rhs
        return self

    @log_call
    def __gt__(self, rhs):
        self.root.append((self, rhs))
        return (self, rhs)

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

        if self.lhs is not None:
            return str(self.lhs) + "." + label
        else:
            return label

class ParseObjSet:

    def __init__(self, objs):
        if len(objs) < 2 or any([not obj.called for obj in objs]):
            raise RuntimeError
        self.objs = list(objs)

    @property
    def root(self):
        return self.objs[-1].root

    def __coerce__(self, other):
        return None

    @log_call
    def __or__(self, rhs):
        self.objs[-1] = (self.objs[-1] | rhs)
        return self

    @log_call
    def __gt__(self, rhs):
        self.root.append((self, rhs))
        return (self, rhs)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        labels = [str(obj) for obj in self.objs]
        return "+".join(labels)

class ParseObjPartial:

    def __init__(self, lhs, name):
        self.lhs = lhs
        self.name = name

    @log_call
    def __call__(self, *args, **kwargs):
        rhs = ParseObj(self.lhs.root, self.name, self.lhs)
        return rhs(*args, **kwargs)

    def __coerce__(self, other):
        return None

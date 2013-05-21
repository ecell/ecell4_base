from logger import log_call


class AnyCallable:

    def __init__(self, root, name, lhs=None):
        self.root = root # a reference to cache

        self.lhs = lhs
        self.name = name
        self.args = []
        self.kwargs = []
        self.key = None
        self.param = None
        self.called = False

    def copy(self):
        retval = AnyCallable(self.root, self.name, self.lhs)
        retval.args = self.args
        retval.kwargs = self.kwargs
        retval.key = self.key
        self.param = self.param
        retval.called = self.called
        return retval

    @log_call
    def __call__(self, *args, **kwargs):
        retval = self.copy()
        retval.args = args
        retval.kwargs = kwargs
        retval.called = True
        return retval

    @log_call
    def __add__(self, rhs):
        return AnyCallablePair(self, rhs)

    def __getitem__(self, key):
        if self.called:
            raise RuntimeError
        retval = self.copy()
        retval.key = key
        return retval

    def __getattr__(self, key):
        return AnyPartial(self, key)

    def __coerce__(self, other):
        return None

    @log_call
    def __or__(self, rhs):
        retval = self.copy()
        retval.param = rhs
        return retval

    @log_call
    def __gt__(self, rhs):
        self.root.append((self, rhs))
        return (self, rhs)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        label = self.name
        if len(self.args) > 0 or len(self.kwargs) > 0:
            attrs = ["%s" for v in self.args]
            attrs += ["%s=%s" % (k, v) for k, v in self.kwargs.items()]
            label += "(%s)" % (",".join(attrs))
        if self.key is not None:
            label += "[%s]" % self.key
        if self.param is not None:
            label += " | %s" % self.param

        if self.lhs is not None:
            return str(self.lhs) + "." + label
        else:
            return label

class AnyCallablePair:

    def __init__(self, lhs, rhs):
        if not lhs.called or not rhs.called:
            raise RuntimeError
        self.lhs = lhs
        self.rhs = rhs

    def __coerce__(self, other):
        return None

    @log_call
    def __or__(self, rhs):
        self.rhs = (self.rhs | rhs)
        return self

    @log_call
    def __gt__(self, rhs):
        self.lhs.root.append((self, rhs)) # XXX: self.lhs does not always have root
        return (self, rhs)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "%s + %s" % (self.lhs, self.rhs)

class AnyPartial:

    def __init__(self, lhs, name):
        self.lhs = lhs
        self.name = name

    @log_call
    def __call__(self, *args, **kwargs):
        rhs = AnyCallable(self.lhs.root, self.name, self.lhs)
        return rhs(*args, **kwargs)

    def __coerce__(self, other):
        return None


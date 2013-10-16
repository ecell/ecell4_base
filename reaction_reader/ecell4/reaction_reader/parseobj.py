import operator
import copy
from logger import log_call


class AnyCallable:
    """AnyCallable must be immutable.
All the members must start with '_'."""

    def __init__(self, root, name):
        self.__root = root # a reference to a callback
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
        self.modification = None

    def set_arguments(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def set_key(self, key):
        self.key = key

    def set_modification(self, rhs):
        self.modification = rhs

    def __str__(self):
        label = self.name

        if self.args is not None or self.kwargs is not None:
            attrs = []
            if self.args is not None:
                attrs += ["%s" % str(v) for v in self.args]
            if self.kwargs is not None:
                attrs += ["%s=%s" % (k, v) for k, v in self.kwargs.items()]
            label += "(%s)" % (",".join(attrs))

        if self.modification is not None:
            label += "^%s" % str(self.modification)
        if self.key is not None:
            label += "[%s]" % str(self.key)
        return label

    def __repr__(self):
        return "<%s.%s: %s>" % (
            self.__class__.__module__, self.__class__.__name__, str(self))

class ExpBase(object):

    def __init__(self, root):
        self.__root = root

    def __repr__(self):
        return "<%s.%s: %s>" % (
            self.__class__.__module__, self.__class__.__name__, str(self))

    def __invert__(self):
        return self.__inv__()

    def __inv__(self):
        retval = InvExp(self.__root, self)
        self.__root.notify_unary_operations(retval)
        return retval

    def __gt__(self, rhs):
        retval = GtExp(self.__root, self, rhs)
        self.__root.notify_comparisons(retval)
        return retval

    def __eq__(self, rhs):
        retval = EqExp(self.__root, self, rhs)
        self.__root.notify_comparisons(retval)
        return retval

    def __ne__(self, rhs):
        retval = NeExp(self.__root, self, rhs)
        self.__root.notify_comparisons(retval)
        return retval

    def __add__(self, rhs):
        retval = AddExp(self.__root, self, rhs)
        return retval

    def __or__(self, rhs):
        retval = OrExp(self.__root, self, rhs)
        self.__root.notify_bitwise_operations(retval)
        return retval

    def __coerce__(self, other):
        return None

class ParseObj(ExpBase):
    """All the members must start with '_'."""

    def __init__(self, root, name, elems=[]):
        ExpBase.__init__(self, root)

        self.__elems = elems + [ParseElem(name)]

    def _elements(self):
        return copy.copy(self.__elems)

    def _size(self):
        return len(self.__elems)

    @log_call
    def __call__(self, *args, **kwargs):
        self.__elems[-1].set_arguments(*args, **kwargs)
        return self

    @log_call
    def __getitem__(self, key):
        self.__elems[-1].set_key(key)
        return self

    @log_call
    def __xor__(self, rhs):
        self.__elems[-1].set_modification(rhs)
        return self

    @log_call
    def __getattr__(self, key):
        if key[0] == "_" and len(key) > 1 and not key[1: ].isdigit():
            raise RuntimeError, (
                "'%s' object has no attribute '%s'"
                    % (self.__class__.__name__, key))
        self.__elems.append(ParseElem(key))
        return self

    def __str__(self):
        labels = [str(elem) for elem in self.__elems]
        return ".".join(labels)

class UnaryExp(ExpBase):

    def __init__(self, root, target, opr=""):
        ExpBase.__init__(self, root)
        self.__target = target
        self.__opr = opr

    @property
    def _target(self):
        return self.__target

    def __str__(self):
        return "%s%s" % (self.__opr, self.__target)

class InvExp(UnaryExp):

    def __init__(self, root, target):
        UnaryExp.__init__(self, root, target, "~")

class AddExp(ExpBase):

    def __init__(self, root, lhs, rhs):
        ExpBase.__init__(self, root)

        self.__elems = []
        self.__append(lhs)
        self.__append(rhs)

    def _elements(self):
        return copy.copy(self.__elems)

    def __append(self, obj):
        if isinstance(obj, AnyCallable):
            self.__elems.append(obj._as_ParseObj())
        elif isinstance(obj, AddExp):
            self.__elems.extend(obj._elements())
        else:
            self.__elems.append(obj)

    def __str__(self):
        return "(%s)" % ("+".join([str(obj) for obj in self.__elems]))

class OrExp(ExpBase):

    def __init__(self, root, lhs, rhs):
        ExpBase.__init__(self, root)

        self.__elems = []
        self.__append(lhs)
        self.__append(rhs)

    def _elements(self):
        return copy.copy(self.__elems)

    def __append(self, obj):
        if isinstance(obj, AnyCallable):
            self.__elems.append(obj._as_ParseObj())
        elif isinstance(obj, OrExp):
            self.__elems.extend(obj._elements())
        else:
            self.__elems.append(obj)

    def __str__(self):
        return "|".join([str(obj) for obj in self.__elems])

class CmpExp(ExpBase):

    def __init__(self, root, lhs, rhs, opr=""):
        ExpBase.__init__(self, root)
        self.__lhs = lhs
        self.__rhs = rhs
        self.__opr = opr

    @property
    def _lhs(self):
        return self.__lhs

    @property
    def _rhs(self):
        return self.__rhs

    def __str__(self):
        return "%s%s%s" % (self.__lhs, self.__opr, self.__rhs)

class GtExp(CmpExp):

    def __init__(self, root, lhs, rhs):
        CmpExp.__init__(self, root, lhs, rhs, ">")

class NeExp(CmpExp):

    def __init__(self, root, lhs, rhs):
        CmpExp.__init__(self, root, lhs, rhs, "<>")

class EqExp(CmpExp):

    def __init__(self, root, lhs, rhs):
        CmpExp.__init__(self, root, lhs, rhs, "==")

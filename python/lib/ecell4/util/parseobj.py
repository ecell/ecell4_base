import operator
import copy
from logger import log_call
import inspect


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

    def __deepcopy__(self, key):
        return AnyCallable(self.__root, self.__name)

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

    def update_arguments(self, *args, **kwargs):
        if self.args is None:
            self.args = args
        else:
            self.args += args

        if self.kwargs is None:
            self.kwargs = kwargs
        else:
            self.kwargs.update(kwargs)

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

    @property
    def _root(self):
        return self.__root

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

    def __init__(self, root, name=None, elems=None):
        ExpBase.__init__(self, root)

        if name is None and elems is None:
            raise RuntimeError, "a name or elements must be given"

        if elems is None:
            self.__elems = []
        else:
            self.__elems = copy.copy(elems)

        if name is not None:
            self.__elems += [ParseElem(name)]

    def _elements(self):
        return copy.copy(self.__elems)

    def _size(self):
        return len(self.__elems)

    def _last(self):
        return self.__elems[-1]

    def _append(self, elem):
        if not isinstance(elem, ParseElem):
            raise RuntimeError, "an invalid unit name was given [%s]" % (
                str(elem))
        self.__elems.append(elem)

    def _extend(self, elems):
        for elem in elems:
            if not isinstance(elem, ParseElem):
                raise RuntimeError, "an invalid unit name was given [%s]" % (
                    str(elem))
        self.__elems.extend(elems)

    def _eval_last(self, *args, **kwargs):
        obj = self.__elems.pop()
        return obj(*args, **kwargs)

    def __deepcopy__(self, memo):
        return ParseObj(self._root, elems=copy.deepcopy(self.__elems))

    @log_call
    def __call__(self, *args, **kwargs):
        if len(args) == 0 and len(kwargs) == 0:
            return self

        retval = copy.deepcopy(self)
        retval._last().update_arguments(*args, **kwargs)
        return retval

    @log_call
    def __getitem__(self, key):
        retval = copy.deepcopy(self)
        retval._last().set_key(key)
        return retval

    @log_call
    def __xor__(self, rhs):
        retval = copy.deepcopy(self)
        retval._last().set_modification(rhs)
        return retval

    # @log_call #XXX: donot wrap
    def __getattr__(self, key):
        if key[0] == "_" and len(key) > 1 and not key[1: ].isdigit():
            raise RuntimeError, (
                "'%s' object has no attribute '%s'"
                    % (self.__class__.__name__, key))

        retval = copy.deepcopy(self)
        calling_frame = inspect.currentframe().f_back

        try:
            obj = eval(key, calling_frame.f_globals, calling_frame.f_locals)
        except NameError:
            retval._append(ParseElem(key))
            return retval

        if isinstance(obj, (AnyCallable, ParseObj)):
            retval._extend(obj._elements())
        elif isinstance(obj, ParseElem):
            retval._append(obj)
        else:
            raise ValueError, (
                "[%s] must be either ParseObj or ParseElem. [%s] given."
                % (key, str(obj)))
        return retval

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

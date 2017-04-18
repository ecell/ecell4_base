import operator
import copy
from .logger import log_call
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

    def __getitem__(self, key):
        return operator.getitem(self._as_ParseObj(), key)

    def __call__(self, *args, **kwargs):
        return self._as_ParseObj()(*args, **kwargs)

    # operators

    def __inv__(self):
        # return operator.inv(self._as_ParseObj())
        retval = InvExp(self.__root, self)
        self.__root.notify_unary_operations(retval)
        return retval

    def __invert__(self):
        return self.__inv__()

    def __pos__(self):
        return operator.pos(self._as_ParseObj())

    def __neg__(self):
        return operator.neg(self._as_ParseObj())

    def __xor__(self, rhs):
        return operator.xor(self._as_ParseObj(), rhs)

    def __gt__(self, rhs):
        return operator.gt(self._as_ParseObj(), rhs)

    def __eq__(self, rhs):
        return operator.eq(self._as_ParseObj(), rhs)

    def __ne__(self, rhs):
        return operator.ne(self._as_ParseObj(), rhs)

    # def __lshift__(self, rhs):
    #     return operator.lshift(self._as_ParseObj(), rhs)

    # def __rlshift__(self, lhs):
    #     return operator.lshift(lhs, self._as_ParseObj())

    # bitwise operations

    def __or__(self, rhs):
        return operator.or_(self._as_ParseObj(), rhs)

    def __ror__(self, lhs):
        return operator.or_(lhs, self._as_ParseObj())

    def __and__(self, rhs):
        return operator.and_(self._as_ParseObj(), rhs)

    def __rand__(self, lhs):
        return operator.and_(lhs, self._as_ParseObj())

    # arithmatic operators

    def __add__(self, rhs):
        return operator.add(self._as_ParseObj(), rhs)

    def __radd__(self, lhs):
        return operator.add(lhs, self._as_ParseObj())

    def __sub__(self, rhs):
        return operator.sub(self._as_ParseObj(), rhs)

    def __rsub__(self, lhs):
        return operator.sub(lhs, self._as_ParseObj())

    def __div__(self, rhs):
        return operator.div(self._as_ParseObj(), rhs)

    def __rdiv__(self, lhs):
        return operator.div(lhs, self._as_ParseObj())

    def __truediv__(self, rhs):
        return operator.truediv(self._as_ParseObj(), rhs)

    def __rtruediv__(self, lhs):
        return operator.truediv(lhs, self._as_ParseObj())

    def __mul__(self, rhs):
        return operator.mul(self._as_ParseObj(), rhs)

    def __rmul__(self, lhs):
        return operator.mul(lhs, self._as_ParseObj())

    def __pow__(self, rhs):
        return operator.pow(self._as_ParseObj(), rhs)

    def __rpow__(self, lhs):
        return operator.pow(lhs, self._as_ParseObj())

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

    def __coerce__(self, other):
        return None

    # operators

    def __inv__(self):
        retval = InvExp(self.__root, self)
        self.__root.notify_unary_operations(retval)
        return retval

    def __invert__(self):
        return self.__inv__()

    def __pos__(self):
        retval = PosExp(self.__root, self)
        self.__root.notify_unary_operations(retval)
        return retval

    def __neg__(self):
        retval = NegExp(self.__root, self)
        self.__root.notify_unary_operations(retval)
        return retval

    # def __xor__(self, rhs):

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

    # def __lshift__(self, rhs):
    #     retval = LshiftExp(self.__root, self, rhs)
    #     return retval

    # def __rlshift__(self, lhs):
    #     retval = LshiftExp(self.__root, lhs, self)
    #     return retval

    # bitwise operations

    def __or__(self, rhs):
        retval = OrExp(self.__root, self, rhs)
        self.__root.notify_bitwise_operations(retval)
        return retval

    def __ror__(self, lhs):
        retval = OrExp(self.__root, lhs, self)
        self.__root.notify_bitwise_operations(retval)
        return retval

    def __and__(self, rhs):
        retval = AndExp(self.__root, self, rhs)
        self.__root.notify_bitwise_operations(retval)
        return retval

    def __rand__(self, lhs):
        retval = AndExp(self.__root, lhs, self)
        self.__root.notify_bitwise_operations(retval)
        return retval

    # operators

    def __add__(self, rhs):
        retval = AddExp(self.__root, self, rhs)
        return retval

    def __radd__(self, lhs):
        retval = AddExp(self.__root, lhs, self)
        return retval

    def __sub__(self, rhs):
        retval = SubExp(self.__root, self, rhs)
        return retval

    def __rsub__(self, lhs):
        retval = SubExp(self.__root, lhs, self)
        return retval

    def __div__(self, rhs):
        retval = DivExp(self.__root, self, rhs)
        return retval

    def __rdiv__(self, lhs):
        retval = DivExp(self.__root, lhs, self)
        return retval

    def __truediv__(self, rhs):
        retval = DivExp(self.__root, self, rhs)
        return retval

    def __rtruediv__(self, lhs):
        retval = DivExp(self.__root, lhs, self)
        return retval

    def __mul__(self, rhs):
        retval = MulExp(self.__root, self, rhs)
        return retval

    def __rmul__(self, lhs):
        retval = MulExp(self.__root, lhs, self)
        return retval

    def __pow__(self, rhs):
        retval = PowExp(self.__root, self, rhs)
        return retval

    def __rpow__(self, lhs):
        retval = PowExp(self.__root, lhs, self)
        return retval

class ParseObj(ExpBase):
    """All the members must start with '_'."""

    def __init__(self, root, name=None, elems=None):
        ExpBase.__init__(self, root)

        if name is None and elems is None:
            raise RuntimeError("a name or elements must be given")

        if elems is None:
            self._elems = []
        else:
            self._elems = copy.copy(elems)

        if name is not None:
            self._elems += [ParseElem(name)]

    def _elements(self):
        return copy.copy(self._elems)

    def _size(self):
        return len(self._elems)

    def _last(self):
        return self._elems[-1]

    def _append(self, elem):
        if not isinstance(elem, ParseElem):
            raise RuntimeError("an invalid unit name was given [%s]" % (
                str(elem)))
        self._elems.append(elem)

    def _extend(self, elems):
        for elem in elems:
            if not isinstance(elem, ParseElem):
                raise RuntimeError("an invalid unit name was given [%s]" % (
                    str(elem)))
        self._elems.extend(elems)

    def _eval_last(self, *args, **kwargs):
        obj = self._elems.pop()
        return obj(*args, **kwargs)

    def __deepcopy__(self, memo):
        return ParseObj(self._root, elems=copy.deepcopy(self._elems))

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
            raise RuntimeError(
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
            raise ValueError(
                "[%s] must be either ParseObj or ParseElem. [%s] given."
                % (key, str(obj)))
        return retval

    def __str__(self):
        labels = [str(elem) for elem in self._elems]
        return ".".join(labels)

class UnaryExp(ExpBase):

    def __init__(self, root, target, opr=""):
        ExpBase.__init__(self, root)
        self._elems = [target]
        self.__opr = opr

    def _elements(self):
        return copy.copy(self.__elements)

    @property
    def _target(self):
        return self._elems[0]

    def __str__(self):
        return "%s%s" % (self.__opr, self._elems[0])

class InvExp(UnaryExp):

    def __init__(self, root, target):
        UnaryExp.__init__(self, root, target, "~")

    def __deepcopy__(self, memo):
        return InvExp(self._root, copy.deepcopy(self._target))

class PosExp(UnaryExp):

    def __init__(self, root, target):
        UnaryExp.__init__(self, root, target, "+")

    def __deepcopy__(self, memo):
        return PosExp(self._root, copy.deepcopy(self._target))

class NegExp(UnaryExp):

    def __init__(self, root, target):
        UnaryExp.__init__(self, root, target, "-")

    def __deepcopy__(self, memo):
        return NegExp(self._root, copy.deepcopy(self._target))

class AddExp(ExpBase):

    def __init__(self, root, lhs, rhs):
        ExpBase.__init__(self, root)

        self._elems = []
        self.__append(lhs)
        self.__append(rhs)

    def _elements(self):
        return copy.copy(self._elems)

    def __append(self, obj):
        if isinstance(obj, AnyCallable):
            self._elems.append(obj._as_ParseObj())
        elif isinstance(obj, AddExp):
            self._elems.extend(obj._elements())
        else:
            self._elems.append(obj)

    def __str__(self):
        return "(%s)" % ("+".join([str(obj) for obj in self._elems]))

    def __deepcopy__(self, memo):
        retval = AddExp(self._root, None, None)
        retval._elems = copy.deepcopy(self._elems)
        return retval

class SubExp(ExpBase):

    def __init__(self, root, lhs, rhs):
        ExpBase.__init__(self, root)

        self._elems = []
        self.__append(lhs)
        self.__append(rhs)

    def _elements(self):
        return copy.copy(self._elems)

    def __append(self, obj):
        if isinstance(obj, AnyCallable):
            self._elems.append(obj._as_ParseObj())
        elif len(self._elems) > 0 and isinstance(obj, SubExp):
            self._elems.extend(obj._elements())
        else:
            self._elems.append(obj)

    def __str__(self):
        return "(%s)" % ("-".join([str(obj) for obj in self._elems]))

    def __deepcopy__(self, memo):
        retval = SubExp(self._root, None, None)
        retval._elems = copy.deepcopy(self._elems)
        return retval

class DivExp(ExpBase):

    def __init__(self, root, lhs, rhs):
        ExpBase.__init__(self, root)

        self._elems = []
        self.__append(lhs)
        self.__append(rhs)

    def _elements(self):
        return copy.copy(self._elems)

    def __append(self, obj):
        if isinstance(obj, AnyCallable):
            self._elems.append(obj._as_ParseObj())
        elif len(self._elems) > 0 and isinstance(obj, DivExp):
            self._elems.extend(obj._elements())
        else:
            self._elems.append(obj)

    def __str__(self):
        return "(%s)" % ("/".join([str(obj) for obj in self._elems]))

    def __deepcopy__(self, memo):
        retval = DivExp(self._root, None, None)
        retval._elems = copy.deepcopy(self._elems)
        return retval

class MulExp(ExpBase):

    def __init__(self, root, lhs, rhs):
        ExpBase.__init__(self, root)

        self._elems = []
        self.__append(lhs)
        self.__append(rhs)

    def _elements(self):
        return copy.copy(self._elems)

    def __append(self, obj):
        if isinstance(obj, AnyCallable):
            self._elems.append(obj._as_ParseObj())
        elif isinstance(obj, MulExp):
            self._elems.extend(obj._elements())
        else:
            self._elems.append(obj)

    def __str__(self):
        return "(%s)" % ("*".join([str(obj) for obj in self._elems]))

    def __deepcopy__(self, memo):
        retval = MulExp(self._root, None, None)
        retval._elems = copy.deepcopy(self._elems)
        return retval

# class LshiftExp(ExpBase):
# 
#     def __init__(self, root, lhs, rhs):
#         ExpBase.__init__(self, root)
# 
#         self._elems = []
#         self.__append(lhs)
#         self.__append(rhs)
# 
#     def _elements(self):
#         return copy.copy(self._elems)
# 
#     def __append(self, obj):
#         if isinstance(obj, AnyCallable):
#             self._elems.append(obj._as_ParseObj())
#         elif isinstance(obj, LshiftExp):
#             self._elems.extend(obj._elements())
#         else:
#             self._elems.append(obj)
# 
#     def __str__(self):
#         return "(%s)" % ("<<".join([str(obj) for obj in self._elems]))

class PowExp(ExpBase):

    def __init__(self, root, lhs, rhs):
        ExpBase.__init__(self, root)
        self._elems = [lhs, rhs]

    @property
    def _lhs(self):
        return self._elems[0]

    @property
    def _rhs(self):
        return self._elems[1]

    def _elements(self):
        return copy.copy(self._elems)

    def __str__(self):
        return "pow(%s,%s)" % (self._elems[0], self._elems[1])
        # return "(%s**%s)" % (self._elems[0], self._elems[1])

    def __deepcopy__(self, memo):
        retval = PowExp(self._root, None, None)
        retval._elems = copy.deepcopy(self._elems)
        return retval

class OrExp(ExpBase):

    def __init__(self, root, lhs, rhs):
        ExpBase.__init__(self, root)

        self._elems = []
        self.__append(lhs)
        self.__append(rhs)

    def _elements(self):
        return copy.copy(self._elems)

    def __append(self, obj):
        if isinstance(obj, AnyCallable):
            self._elems.append(obj._as_ParseObj())
        elif isinstance(obj, OrExp):
            self._elems.extend(obj._elements())
        else:
            self._elems.append(obj)

    def __str__(self):
        return "(%s)" % ("|".join([str(obj) for obj in self._elems]))

    def __deepcopy__(self, memo):
        retval = OrExp(self._root, None, None)
        retval._elems = copy.deepcopy(self._elems)
        return retval

class AndExp(ExpBase):

    def __init__(self, root, lhs, rhs):
        ExpBase.__init__(self, root)

        self._elems = []
        self.__append(lhs)
        self.__append(rhs)

    def _elements(self):
        return copy.copy(self._elems)

    def __append(self, obj):
        if isinstance(obj, AnyCallable):
            self._elems.append(obj._as_ParseObj())
        elif isinstance(obj, AndExp):
            self._elems.extend(obj._elements())
        else:
            self._elems.append(obj)

    def __str__(self):
        return "(%s)" % ("&".join([str(obj) for obj in self._elems]))

    def __deepcopy__(self, memo):
        retval = AndExp(self._root, None, None)
        retval._elems = copy.deepcopy(self._elems)
        return retval

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

    def __deepcopy__(self, memo):
        return GtExp(self._root, copy.deepcopy(self._lhs), copy.deepcopy(self._rhs))

class NeExp(CmpExp):

    def __init__(self, root, lhs, rhs):
        CmpExp.__init__(self, root, lhs, rhs, "<>")

    def __deepcopy__(self, memo):
        return NeExp(self._root, copy.deepcopy(self._lhs), copy.deepcopy(self._rhs))

class EqExp(CmpExp):

    def __init__(self, root, lhs, rhs):
        CmpExp.__init__(self, root, lhs, rhs, "==")

    def __deepcopy__(self, memo):
        return EqExp(self._root, copy.deepcopy(self._lhs), copy.deepcopy(self._rhs))

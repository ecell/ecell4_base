import textwrap
from logging import getLogger
_log = getLogger(__name__)

from ..util.parseobj import AnyCallable, ExpBase

import pint
from pint.quantity import _Quantity
from pint.unit import _Unit
from pint.errors import UndefinedUnitError

__all__ = [
    'getUnitRegistry', '_Quantity', '_Unit', 'check_dimensionality', 'wrap_quantity',
    'get_application_registry']

def wrapped_binary_operator(op1, op2):
    def wrapped(self, other):
        if isinstance(other, ExpBase):
            return op2(other, self)
        elif isinstance(other, AnyCallable):
            return op2(other._as_ParseObj(), self)
        return op1(self, other)
    return wrapped

def wrap_quantity(cls):
    cls.__add__ = wrapped_binary_operator(cls.__add__, ExpBase.__radd__)
    cls.__mul__ = wrapped_binary_operator(cls.__mul__, ExpBase.__rmul__)
    cls.__div__ = wrapped_binary_operator(cls.__div__, ExpBase.__rdiv__)
    cls.__truediv__ = wrapped_binary_operator(cls.__truediv__, ExpBase.__rtruediv__)
    cls.__pow__ = wrapped_binary_operator(cls.__pow__, ExpBase.__rpow__)
    return cls

def getUnitRegistry(length="meter", time="second", substance="item", volume=None, other=()):
    """Return a pint.UnitRegistry made compatible with ecell4.

    Parameters
    ----------
    length : str, optional
        A default unit for '[length]'. 'meter' is its default.
    time : str, optional
        A default unit for '[time]'. 'second' is its default.
    substance : str, optional
        A default unit for '[substance]' (the number of molecules). 'item' is its default.
    volume : str, optional
        A default unit for '[volume]'. Its default is None, thus '[length]**3'.
    other : tuple, optional
        A list of user-defined default units other than the above.

    Returns
    -------
    ureg : pint.UnitRegistry

    """
    ureg = pint.UnitRegistry()
    ureg.define('item = mole / (avogadro_number * 1 mole)')

    try:
        pint.molar
    # except UndefinedUnitError:
    except AttributeError:
        # https://github.com/hgrecco/pint/blob/master/pint/default_en.txt#L75-L77
        ureg.define('[concentration] = [substance] / [volume]')
        ureg.define('molar = mol / (1e-3 * m ** 3) = M')

    base_units = [unit for unit in (length, time, substance, volume) if unit is not None]
    base_units.extend(other)
    s = ureg.System.from_lines(
        ["@system local using international"] + base_units,
        ureg.get_base_units)
    ureg.default_system = 'local'

    wrap_quantity(ureg.Quantity)

    pint.set_application_registry(ureg)  # for pickling
    return ureg

def check_dimensionality(q, dim):
    """Return whether the quantity has the dimensionality given.

    Parameters
    ----------
    q : pint._Quantity
    dim : pint.util.UnitsContainer

    Returns
    -------
    value : bool

    """
    return (q._REGISTRY.get_dimensionality(dim) == q.dimensionality)

def get_application_registry():
    """Just return `pint._APP_REGISTRY`."""
    return pint._APP_REGISTRY


if __name__ == '__main__':
    ureg = getUnitRegistry(length='micrometer', volume='liter')
    Q_ = ureg.Quantity
    q = Q_(1.0, "ml")
    print(q)
    print(q.to_base_units())
    q = Q_(1.0, "M")
    print(q)
    print(q.to_base_units())

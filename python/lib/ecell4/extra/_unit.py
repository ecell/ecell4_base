import textwrap
from logging import getLogger
_log = getLogger(__name__)

from ..util.parseobj import AnyCallable, ExpBase

import pint
from pint.quantity import _Quantity
from pint.errors import UndefinedUnitError

__all__ = ['getUnitRegistry', '_Quantity']

def wrapped_binary_operator(op1, op2):
    def wrapped(self, other):
        if isinstance(other, ExpBase):
            print('wrapped:', self, other, type(self), type(other))
            return op2(other, self)
        elif isinstance(other, AnyCallable):
            print('wrapped:', self, other, type(self), type(other))
            return op2(other._as_ParseObj(), self)
        return op1(self, other)
    return wrapped

def getUnitRegistry(length="meter", time="second", substance="item", other=()):
    ureg = pint.UnitRegistry()
    ureg.define('item = mole / (avogadro_number * 1 mole)')

    try:
        pint.molar
    # except UndefinedUnitError:
    except AttributeError:
        # https://github.com/hgrecco/pint/blob/master/pint/default_en.txt#L75-L77
        ureg.define('[concentration] = [substance] / [volume]')
        ureg.define('molar = mol / (1e-3 * m ** 3) = M')

    conf = textwrap.dedent("""\
        @system local using international
        {}
        {}
        {}""".format(length, time, substance)).split('\n')
    conf.extend(other)
    s = ureg.System.from_lines(conf, ureg.get_base_units)
    ureg.default_system = 'local'

    ureg.Quantity.__add__ = wrapped_binary_operator(ureg.Quantity.__add__, ExpBase.__radd__)
    ureg.Quantity.__mul__ = wrapped_binary_operator(ureg.Quantity.__mul__, ExpBase.__rmul__)
    ureg.Quantity.__div__ = wrapped_binary_operator(ureg.Quantity.__div__, ExpBase.__rdiv__)
    ureg.Quantity.__truediv__ = wrapped_binary_operator(ureg.Quantity.__truediv__, ExpBase.__rtruediv__)
    ureg.Quantity.__pow__ = wrapped_binary_operator(ureg.Quantity.__pow__, ExpBase.__rpow__)

    return ureg

# def Quantity(value, unit):
#     q = value * ureg.parse_units(unit)
#     _log.debug("{} => {}".format(q, q.to_base_units()))  #XXX: Just for debugging
#     return q.to_base_units().magnitude
# 
# Q_ = Quantity


if __name__ == '__main__':
    ureg = getUnitRegistry(length='micrometer')
    Q_ = ureg.Quantity
    q = Q_(1.0, "ml")
    print(q)
    print(q.to_base_units())

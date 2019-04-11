from ecell4_base.core import *

def setUpEqualities(self):
    self.addTypeEqualityFunc(Quantity_Real, _assertEqualsQuantity(self))
    self.addTypeEqualityFunc(Quantity_Integer, _assertEqualsQuantity(self))

def _assertEqualsQuantity(self):
    def wrapper(x, y, msg=None):
        self.assertEqual(x.magnitude, y.magnitude, msg)
        self.assertEqual(x.units, y.units, msg)
    return wrapper

import unittest
from ecell4_base.core import *
import math

class Real3Test(unittest.TestCase):
    def setUp(self):
        self.addTypeEqualityFunc(Real3, Real3.__eq__)

    def test_constructor(self):
        x = Real3(1, 2, 3)

    def test_tuple_and_list(self):
        x = Real3(1, 2, 3)
        self.assertEqual(tuple(x), (1.0, 2.0, 3.0))
        self.assertEqual(list(x), [1.0, 2.0, 3.0])

    def test_operators(self):
        x = Real3(7, 4, 9)
        y = Real3(1, 2, 3)
        self.assertEqual(x + y, Real3(8, 7, 12))
        self.assertEqual(x - y, Real3(6, 2, 6))
        self.assertEqual(x * 2, Real3(14, 8, 18))
        self.assertEqual(2 * x, Real3(14, 8, 18))
        self.assertEqual(x / 3, Real3(3.5, 2, 4.5))

    def test_abs(self):
        x = Real3(1, 2, 3)
        self.assertEqual(abs(x), x)
        self.assertEqual(abs(Real3(-1, 2, -3)), Real3(1, 2, 3))

    def test_length(self):
        x = Real3(1, 2, 3)
        sq = 1*1 + 2*2 + 3*3
        self.assertEqual(length_sq(x), sq)
        self.assertEqual(length(x), math.sqrt(sq))

    def test_product(self):
        x = Real3(1, 2, 3)
        y = Real3(4, 5, 2)
        self.assertEqual(dot_product(x, y), 1*4 + 2*5 + 3*2)
        self.assertEqual(dot_product(y, x), 1*4 + 2*5 + 3*2)
        self.assertEqual(cross_product(x, y), Real3(-11, 10, -3))
        self.assertEqual(cross_product(x, y), Real3(11, -10, 3))

class Integer3Test(unittest.TestCase):
    def setUp(self):
        self.addTypeEqualityFunc(Integer3, Integer3.__eq__)

    def test_constructor(self):
        x = Integer3(1, 2, 3)

    def test_tuple_and_list(self):
        x = Integer3(1, 2, 3)
        self.assertEqual(tuple(x), (1, 2, 3))
        self.assertEqual(list(x), [1, 2, 3])

    def test_operators(self):
        x = Integer3(6, 0, 2)
        y = Integer3(1, 2, 3)
        self.assertEqual(x + y, Integer3(7, 2, 5))
        self.assertEqual(x - y, Integer3(5, -2, 1))
        self.assertEqual(x * 2, Integer3(12, 0, 4))
        # self.assertEqual(2 * x, Integer3(12, 0, 4))

    def test_abs(self):
        self.assertEqual(abs(Integer3(1, 2, 3)), Integer3(1, 2, 3))
        self.assertEqual(abs(Integer3(-1, 2, -3)), Integer3(1, 2, 3))

    def test_length(self):
        x = Integer3(1, 2, 3)
        sq = 1*1 + 2*2 + 3*3
        self.assertEqual(length_sq(x), sq)
        self.assertEqual(length(x), math.sqrt(sq))

    def test_product(self):
        x = Integer3(1, 2, 3)
        y = Integer3(4, 5, 2)
        self.assertEqual(dot_product(x, y), 1*4 + 2*5 + 3*2)
        self.assertEqual(dot_product(y, x), 1*4 + 2*5 + 3*2)


class QuantityTest(unittest.TestCase):
    def test_quantity_real(self):
        q = Quantity_Real(0.0)
        self.assertEqual(q.magnitude, 0.0)
        self.assertEqual(q.units, "")

        q = Quantity_Real(1.0, "nm")
        self.assertEqual(q.magnitude, 1.0)
        self.assertEqual(q.units, "nm")

        q.magnitude = 2.0
        q.units = "m/s"
        self.assertEqual(q.magnitude, 2.0)
        self.assertEqual(q.units, "m/s")

    def test_quantity_integer(self):
        q = Quantity_Integer(0)
        self.assertEqual(q.magnitude, 0)
        self.assertEqual(q.units, "")

        q = Quantity_Integer(1, "a")
        self.assertEqual(q.magnitude, 1)
        self.assertEqual(q.units, "a")

        q.magnitude = 2
        q.units = "yen"
        self.assertEqual(q.magnitude, 2)
        self.assertEqual(q.units, "yen")

from ecell4.core import *

import unittest


class SpeciesTest(unittest.TestCase):

    def setUp(self):
        pass

    def test1(self):
        sp = Species('A')
        self.assertEqual(sp.serial(), 'A')

    def test2(self):
        sp = Species()
        sp.add_unit(UnitSpecies('B'))
        sp.add_unit(UnitSpecies('C'))
        sp.add_unit(UnitSpecies('A'))
        self.assertEqual(sp.serial(), 'B.C.A')
        self.assertEqual(sp.num_units(), 3)

    def test3(self):
        sp = Species('A')

        sp.set_attribute('foo', 'bar')
        sp.set_attribute('spam', 'ham')
        sp.set_attribute('hoge', 'hage')

        self.assertTrue(sp.has_attribute('spam'))
        self.assertTrue(sp.has_attribute('foo'))
        self.assertTrue(sp.has_attribute('hoge'))
        self.assertFalse(sp.has_attribute('eggs'))

        sp.remove_attribute('spam')
        self.assertFalse(sp.has_attribute('spam'))

        self.assertEqual(sp.get_attribute('foo'), 'bar')
        self.assertEqual(sp.get_attribute('hoge'), 'hage')

        attrs = sp.list_attributes()
        self.assertEqual(len(attrs), 2)
        for key, value in attrs:
            self.assertTrue(key == 'foo' or key == 'hoge')
            self.assertTrue(
                (key == 'foo' and value == 'bar')
                or (key == 'hoge' and value == 'hage'))

    def test4(self):
        sp = Species()
        sp.deserialize('A.B.C')
        self.assertEqual(sp.serial(), 'A.B.C')
        self.assertEqual(sp.num_units(), 3)

        sp.add_unit(UnitSpecies('D'))
        self.assertEqual(sp.serial(), 'A.B.C.D')
        self.assertEqual(sp.num_units(), 4)

        units = sp.units()
        self.assertEqual(len(units), 4)

    def test5(self):
        sp = Species('X(a,b=c^1).Y(d=e^1,f=g)')
        units = sp.units()

        self.assertEqual(sp.num_units(), 2)
        self.assertEqual(len(units), 2)

        self.assertEqual(units[0].name(), 'X')
        self.assertEqual(units[1].name(), 'Y')

        units[1].add_site('h', 'i', '')

    def test6(self):
        sp = Species(" A   . B . C.D")
        units = sp.units()
        self.assertEqual(len(units), 4)
        self.assertEqual(units[0].name(), "A")
        self.assertEqual(units[1].name(), "B")


if __name__ == '__main__':
    unittest.main()

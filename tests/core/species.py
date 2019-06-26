import unittest
import copy
from ecell4_base.core import *

class UnitSpeciesTest(unittest.TestCase):
    def test_constructor(self):
        usp = UnitSpecies()
        usp = UnitSpecies("A")

    def test_getter(self):
        usp = UnitSpecies("A")
        self.assertEqual(usp.name(), "A")
        self.assertEqual(usp.serial(), "A")


class SpeciesTest(unittest.TestCase):
    def test_serial(self):
        sp = Species('A')
        self.assertEqual(sp.serial(), 'A')

    def test_units(self):
        sp = Species()
        sp.add_unit(UnitSpecies('B'))
        sp.add_unit(UnitSpecies('C'))
        sp.add_unit(UnitSpecies('A'))
        self.assertEqual(sp.serial(), 'B.C.A')
        self.assertEqual(len(sp.units()), 3)

        sp = Species('A.B.C')
        self.assertEqual(sp.serial(), 'A.B.C')
        self.assertEqual(len(sp.units()), 3)

        sp.add_unit(UnitSpecies('D'))
        self.assertEqual(sp.serial(), 'A.B.C.D')
        self.assertEqual(len(sp.units()), 4)

        units = sp.units()
        self.assertEqual(len(units), 4)

        sp = Species('X(a,b=c^1).Y(d=e^1,f=g)')
        units = sp.units()

        self.assertEqual(len(sp.units()), 2)
        self.assertEqual(len(units), 2)

        self.assertEqual(units[0].name(), 'X')
        self.assertEqual(units[1].name(), 'Y')

        units[1].add_site('h', 'i', '')

        sp = Species(" A   . B . C.D")
        units = sp.units()
        self.assertEqual(len(units), 4)
        self.assertEqual(units[0].name(), "A")
        self.assertEqual(units[1].name(), "B")

    def test_attributes(self):
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

    def test_attribute_types(self):
        sp = Species('A')

        sp.set_attribute('key', 'value')
        self.assertTrue(isinstance(sp.get_attribute('key'), str))
        sp.set_attribute('key', True)
        self.assertTrue(isinstance(sp.get_attribute('key'), bool))
        sp.set_attribute('key', Quantity_Integer(2, 'dimensionless'))
        self.assertTrue(isinstance(sp.get_attribute('key'), Quantity_Integer))
        sp.set_attribute('key', Quantity_Real(1.5, 'm**2/s'))
        self.assertTrue(isinstance(sp.get_attribute('key'), Quantity_Real))

        sp.set_attribute('key', 2)
        self.assertTrue(isinstance(sp.get_attribute('key'), Quantity_Integer))
        sp.set_attribute('key', 2.0)
        self.assertTrue(isinstance(sp.get_attribute('key'), Quantity_Real))

    def test_operators(self):
        self.assertTrue(Species('A') == Species('A'))
        self.assertFalse(Species('A') == Species('B'))
        self.assertTrue(Species('A') != Species('B'))

        self.assertTrue(Species('A') < Species('B'))
        self.assertTrue(Species('B') > Species('A'))
        self.assertFalse(Species('A') < Species('A'))

        # self.assertTrue(Species('A') <= Species('A'))   # Not implemented yet

    def test_count_species_matches(self):
        sp = Species("A")
        self.assertEqual(count_species_matches(sp, Species("A")), 1)
        self.assertEqual(count_species_matches(sp, Species("A.A")), 2)

        sp = Species("A.B")
        self.assertEqual(count_species_matches(sp, Species("A.B")), 1)
        self.assertEqual(count_species_matches(sp, Species("B.A")), 1)

        sp = Species("A(p=u^_)")
        self.assertEqual(count_species_matches(sp, Species("A(p=u^1).B(b^1)")), 1)

    def test_pickling(self):
        sp = Species("A")
        sp.set_attribute("key1", "value1")
        sp.set_attribute("key2", Quantity_Real(2.0, "units"))
        sp.set_attribute("key3", True)

        self.assertTrue(sp.has_attribute("key1"))
        self.assertEqual(sp.get_attribute("key1"), "value1")
        self.assertTrue(sp.has_attribute("key2"))
        # self.assertEqual(sp.get_attribute("key2").magnitude, Quantity_Real(2.0, "units").magnitude)
        # self.assertEqual(sp.get_attribute("key2").units, Quantity_Real(2.0, "units").units)
        self.assertEqual(sp.get_attribute("key2"), Quantity_Real(2.0, "units"))
        self.assertTrue(sp.has_attribute("key3"))
        self.assertEqual(sp.get_attribute("key3"), True)

        another = copy.copy(sp)
        self.assertTrue(sp is not another)
        self.assertEqual(sp, another)

        self.assertTrue(another.has_attribute("key1"))
        self.assertEqual(another.get_attribute("key1"), "value1")
        self.assertTrue(another.has_attribute("key2"))
        # self.assertEqual(another.get_attribute("key2").magnitude, Quantity_Real(2.0, "units").magnitude)
        # self.assertEqual(another.get_attribute("key2").units, Quantity_Real(2.0, "units").units)
        self.assertEqual(another.get_attribute("key2"), Quantity_Real(2.0, "units"))
        self.assertTrue(another.has_attribute("key3"))
        self.assertEqual(another.get_attribute("key3"), True)

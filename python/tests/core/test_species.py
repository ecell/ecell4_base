from ecell4.core import *

import unittest


class SpeciesTest(unittest.TestCase):

    def setUp(self):
        pass

    def test1(self):
        sp = Species("A")
        self.assertEqual(sp.serial(), "A")

    def test2(self):
        sp = Species()
        sp.add_unit(UnitSpecies("B"))
        sp.add_unit(UnitSpecies("C"))
        sp.add_unit(UnitSpecies("A"))
        self.assertEqual(sp.serial(), "A.B.C")


if __name__ == '__main__':
    unittest.main()

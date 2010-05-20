#!/usr/bin/env python
import _gfrd
import unittest

class ModelTestCase(unittest.TestCase):
    def setUp(self):
        self.m = _gfrd.Model()

    def tearDown(self):
        pass

    def test_network_rules(self):
        self.assertTrue(isinstance(self.m.network_rules, _gfrd.NetworkRules))

    def test_add_species_type(self):
        s1 = _gfrd.SpeciesType()
        try:
            s1.id
        except _gfrd.IllegalState:
            self.assertTrue(True)
        except:
            self.fail()

    def test_get_species_type_by_id(self):
        s1 = _gfrd.SpeciesType()
        s2 = _gfrd.SpeciesType()
        self.m.add_species_type(s1)
        self.m.add_species_type(s2)
        self.assertNotEqual(s1, s2)
        self.assertNotEqual(s1.id, s2.id)
        self.assertTrue(self.m.get_species_type_by_id(s1.id), s1)
        self.assertTrue(self.m.get_species_type_by_id(s2.id), s2)
        self.assertTrue(self.m.get_species_type_by_id(s1.id).id, s1.id)
        self.assertTrue(self.m.get_species_type_by_id(s2.id).id, s2.id)

if __name__ == "__main__":
    unittest.main()


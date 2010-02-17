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

    def test_new_species_type(self):
        s1 = self.m.new_species_type()
        self.assertNotEqual(None, s1)
        self.assertTrue(isinstance(s1, _gfrd.SpeciesType))
        s2 = self.m.new_species_type()
        self.assertTrue(isinstance(s2, _gfrd.SpeciesType))
        self.assertNotEqual(None, s2)
        self.assertNotEqual(s1, s2)

    def test_get_species_type_by_id(self):
        s1 = self.m.new_species_type()
        s2 = self.m.new_species_type()
        self.assertNotEqual(s1, s2)
        self.assertNotEqual(s1.id, s2.id)
        self.assertTrue(self.m.get_species_type_by_id(s1.id), s1)
        self.assertTrue(self.m.get_species_type_by_id(s2.id), s2)
        self.assertTrue(self.m.get_species_type_by_id(s1.id).id, s1.id)
        self.assertTrue(self.m.get_species_type_by_id(s2.id).id, s2.id)

if __name__ == "__main__":
    unittest.main()


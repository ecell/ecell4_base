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

if __name__ == "__main__":
    unittest.main()


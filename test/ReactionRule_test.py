#!/usr/bin/env python
import _gfrd
import unittest

class ReactionRuleTestCase(unittest.TestCase):
    def setUp(self):
        self.m = _gfrd.Model()
        self.s1 = self.m.new_species_type()
        self.s2 = self.m.new_species_type()

    def tearDown(self):
        pass

    def test_instantiation(self):
        s1, s2 = self.s1, self.s2
        self.assertTrue(isinstance(
            _gfrd.ReactionRule([s1], []), _gfrd.ReactionRule))
        self.assertTrue(isinstance(
            _gfrd.ReactionRule([s1, s2], []), _gfrd.ReactionRule))
        self.assertTrue(isinstance(
            _gfrd.ReactionRule([s1], [s1]), _gfrd.ReactionRule))
        self.assertTrue(isinstance(
            _gfrd.ReactionRule([s1, s1], [s1]), _gfrd.ReactionRule))
        self.assertTrue(isinstance(
            _gfrd.ReactionRule([s1], [s2]), _gfrd.ReactionRule))
        self.assertTrue(isinstance(
            _gfrd.ReactionRule([s1, s1], [s2]), _gfrd.ReactionRule))
        self.assertTrue(isinstance(
            _gfrd.ReactionRule([s1, s1], [s1, s1]), _gfrd.ReactionRule))
        self.assertTrue(isinstance(
            _gfrd.ReactionRule([s1, s1], [s2, s2]), _gfrd.ReactionRule))
        self.assertRaises(TypeError, lambda:
            self.assertTrue(isinstance(
                 _gfrd.ReactionRule([], []), _gfrd.ReactionRule)))
        self.assertRaises(TypeError, lambda:
            self.assertTrue(isinstance(
                _gfrd.ReactionRule([], [s1]), _gfrd.ReactionRule)))
        self.assertRaises(TypeError, lambda:
            self.assertTrue(isinstance(
                _gfrd.ReactionRule([], [s1, s2]), _gfrd.ReactionRule)))

    def test_comparison(self):
        s1, s2 = self.s1, self.s2
        self.assertEqual(
            _gfrd.ReactionRule([s1], []),
            _gfrd.ReactionRule([s1], []))
        self.assertEqual(
            _gfrd.ReactionRule([s1], [s1]),
            _gfrd.ReactionRule([s1], [s1]))
        self.assertEqual(
            _gfrd.ReactionRule([s1], [s2]),
            _gfrd.ReactionRule([s1], [s2]))
        self.assertEqual(
            _gfrd.ReactionRule([s1], [s1, s2]),
            _gfrd.ReactionRule([s1], [s1, s2]))
        self.assertEqual(
            _gfrd.ReactionRule([s1], [s2, s1]),
            _gfrd.ReactionRule([s1], [s2, s1]))
        self.assertEqual(
            _gfrd.ReactionRule([s1], [s1, s2]),
            _gfrd.ReactionRule([s1], [s2, s1]))
        self.assertEqual(
            _gfrd.ReactionRule([s1], [s2, s1]),
            _gfrd.ReactionRule([s1], [s1, s2]))
        self.assertNotEqual(
            _gfrd.ReactionRule([s1], []),
            _gfrd.ReactionRule([s2], []))
        self.assertNotEqual(
            _gfrd.ReactionRule([s1], [s1]),
            _gfrd.ReactionRule([s2], [s1]))
        self.assertNotEqual(
            _gfrd.ReactionRule([s1], [s2]),
            _gfrd.ReactionRule([s2], [s2]))
        self.assertNotEqual(
            _gfrd.ReactionRule([s1], [s1, s2]),
            _gfrd.ReactionRule([s2], [s1, s2]))
        self.assertNotEqual(
            _gfrd.ReactionRule([s1], [s2, s1]),
            _gfrd.ReactionRule([s2], [s2, s1]))
        self.assertNotEqual(
            _gfrd.ReactionRule([s1], [s1, s2]),
            _gfrd.ReactionRule([s2], [s2, s1]))
        self.assertNotEqual(
            _gfrd.ReactionRule([s1], [s2, s1]),
            _gfrd.ReactionRule([s2], [s1, s2]))

    def test_get_attribute(self):
        rr = _gfrd.ReactionRule([self.s1], [])
        rr['k'] = '0.0'
        self.assertEqual('0.0', rr['k'])
        rr['k'] = '0.5'
        self.assertEqual('0.5', rr['k'])
        rr['name'] = 'R1'
        self.assertEqual('R1', rr['name'])

    def test_get_reactants_and_get_products(self):
        s1, s2 = self.s1, self.s2
        for reactants in [(s1, ), (s2, ), (s1, s2),]:
            for products in [(), (s1, ), (s2, ), (s1, s2),]:
                r = _gfrd.ReactionRule(reactants, products)
                self.assertEqual(reactants, tuple(reactants))
                self.assertEqual(products, tuple(products))


if __name__ == "__main__":
    unittest.main()


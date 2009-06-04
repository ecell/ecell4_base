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
            _gfrd.ReactionRule([s1], [], .0), _gfrd.ReactionRule))
        self.assertTrue(isinstance(
            _gfrd.ReactionRule([s1, s2], [], .0), _gfrd.ReactionRule))
        self.assertTrue(isinstance(
            _gfrd.ReactionRule([s1], [s1], .0), _gfrd.ReactionRule))
        self.assertTrue(isinstance(
            _gfrd.ReactionRule([s1, s1], [s1], .0), _gfrd.ReactionRule))
        self.assertTrue(isinstance(
            _gfrd.ReactionRule([s1], [s2], .0), _gfrd.ReactionRule))
        self.assertTrue(isinstance(
            _gfrd.ReactionRule([s1, s1], [s2], .0), _gfrd.ReactionRule))
        self.assertTrue(isinstance(
            _gfrd.ReactionRule([s1, s1], [s1, s1], .0), _gfrd.ReactionRule))
        self.assertTrue(isinstance(
            _gfrd.ReactionRule([s1, s1], [s2, s2], .0), _gfrd.ReactionRule))
        self.assertRaises(TypeError, lambda:
            self.assertTrue(isinstance(
                 _gfrd.ReactionRule([], [], .0), _gfrd.ReactionRule)))
        self.assertRaises(TypeError, lambda:
            self.assertTrue(isinstance(
                _gfrd.ReactionRule([], [s1], .0), _gfrd.ReactionRule)))
        self.assertRaises(TypeError, lambda:
            self.assertTrue(isinstance(
                _gfrd.ReactionRule([], [s1, s2], .0), _gfrd.ReactionRule)))

    def tesT_comparison(self):
        self.assertEqual(
            _gfrd.ReactionRule([s1], [], .0),
            _gfrd.ReactionRule([s1], [], .0))
        self.assertEqual(
            _gfrd.ReactionRule([s1], [s1], .0),
            _gfrd.ReactionRule([s1], [s1], .0))
        self.assertEqual(
            _gfrd.ReactionRule([s1], [s2], .0),
            _gfrd.ReactionRule([s1], [s2], .0))
        self.assertEqual(
            _gfrd.ReactionRule([s1], [s1, s2], .0),
            _gfrd.ReactionRule([s1], [s1, s2], .0))
        self.assertEqual(
            _gfrd.ReactionRule([s1], [s2, s1], .0),
            _gfrd.ReactionRule([s1], [s2, s1], .0))
        self.assertEqual(
            _gfrd.ReactionRule([s1], [s1, s2], .0),
            _gfrd.ReactionRule([s1], [s2, s1], .0))
        self.assertEqual(
            _gfrd.ReactionRule([s1], [s2, s1], .0),
            _gfrd.ReactionRule([s1], [s1, s2], .0))
        self.assertNotEqual(
            _gfrd.ReactionRule([s1], [], .0),
            _gfrd.ReactionRule([s2], [], .0))
        self.assertEqual(
            _gfrd.ReactionRule([s1], [s1], .0),
            _gfrd.ReactionRule([s2], [s1], .0))
        self.assertEqual(
            _gfrd.ReactionRule([s1], [s2], .0),
            _gfrd.ReactionRule([s2], [s2], .0))
        self.assertEqual(
            _gfrd.ReactionRule([s1], [s1, s2], .0),
            _gfrd.ReactionRule([s2], [s1, s2], .0))
        self.assertEqual(
            _gfrd.ReactionRule([s1], [s2, s1], .0),
            _gfrd.ReactionRule([s2], [s2, s1], .0))
        self.assertEqual(
            _gfrd.ReactionRule([s1], [s1, s2], .0),
            _gfrd.ReactionRule([s2], [s2, s1], .0))
        self.assertEqual(
            _gfrd.ReactionRule([s1], [s2, s1], .0),
            _gfrd.ReactionRule([s2], [s1, s2], .0))

    def test_get_k(self):
        self.assertAlmostEqual(.0, _gfrd.ReactionRule([self.s1], [], .0).k)
        self.assertAlmostEqual(.5, _gfrd.ReactionRule([self.s1], [], .5).k)
        self.assertAlmostEqual(1., _gfrd.ReactionRule([self.s1], [], 1.).k)

    def test_get_reactants_and_get_products(self):
        s1, s2 = self.s1, self.s2
        for reactants in [(s1, ), (s2, ), (s1, s2),]:
            for products in [(), (s1, ), (s2, ), (s1, s2),]:
                r = _gfrd.ReactionRule(reactants, products, .0)
                self.assertEqual(reactants, tuple(reactants))
                self.assertEqual(products, tuple(products))


if __name__ == "__main__":
    unittest.main()


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
        self.assertRaises(TypeError, lambda:
            self.assertTrue(isinstance(
                 _gfrd.ReactionRule([], [], .0), _gfrd.ReactionRule)))
        self.assertRaises(TypeError, lambda:
            self.assertTrue(isinstance(
                _gfrd.ReactionRule([], [s1], .0), _gfrd.ReactionRule)))
        self.assertRaises(TypeError, lambda:
            self.assertTrue(isinstance(
                _gfrd.ReactionRule([], [s1, s2], .0), _gfrd.ReactionRule)))

    def test_get_k(self):
        self.assertAlmostEqual(.0, _gfrd.ReactionRule([self.s1], [], .0).k)
        self.assertAlmostEqual(.5, _gfrd.ReactionRule([self.s1], [], .5).k)
        self.assertAlmostEqual(1., _gfrd.ReactionRule([self.s1], [], 1.).k)
        

if __name__ == "__main__":
    unittest.main()


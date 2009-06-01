import _gfrd
import unittest

class NetworkRulesTestCase(unittest.TestCase):
    def setUp(self):
        self.m = _gfrd.Model()
        self.s1 = self.m.new_species_type()
        self.s2 = self.m.new_species_type()

    def tearDown(self):
        pass

    def test_add_reaction_rule(self):
        self.m.network_rules.add_reaction_rule(
            _gfrd.ReactionRule([self.s1], [self.s1, self.s2], .2))
        self.assertTrue(True)

        self.m.network_rules.add_reaction_rule(
            _gfrd.ReactionRule([self.s2], [self.s1, self.s2], .2))
        self.assertTrue(True)

        self.assertRaises(_gfrd.AlreadyExists,
                lambda: self.m.network_rules.add_reaction_rule(
                    _gfrd.ReactionRule([self.s1], [self.s1, self.s2], .2)))

        self.assertRaises(_gfrd.AlreadyExists,
                lambda: self.m.network_rules.add_reaction_rule(
                    _gfrd.ReactionRule([self.s2], [self.s1, self.s2], .2)))

    def test_query_reaction_rule(self):
        r1 = _gfrd.ReactionRule([self.s1], [self.s1, self.s2], .2)
        self.m.network_rules.add_reaction_rule(r1)
        a = self.m.network_rules.query_reaction_rule(self.s1)
        self.assertTrue(iter(a) != None)
        a = list(a)
        self.assertEqual(1, len(a))
        self.assertTrue(r1 in a)

        r2 = _gfrd.ReactionRule([self.s1], [self.s1], .2)
        self.m.network_rules.add_reaction_rule(r2)
        a = self.m.network_rules.query_reaction_rule(self.s1)
        self.assertTrue(iter(a) != None)
        a = list(a)
        self.assertEqual(2, len(a))
        self.assertTrue(r1 in a)
        self.assertTrue(r2 in a)

if __name__ == "__main__":
    unittest.main()


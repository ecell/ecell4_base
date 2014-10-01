#!/usr/bin/env python
import _gfrd
import unittest

class NetworkRulesTestCase(unittest.TestCase):
    def setUp(self):
        self.m = _gfrd.Model()
        self.s1 = _gfrd.SpeciesType()
        self.m.add_species_type(self.s1)
        self.s2 = _gfrd.SpeciesType()
        self.m.add_species_type(self.s2)

    def tearDown(self):
        pass

    def test_add_reaction_rule(self):
        self.m.network_rules.add_reaction_rule(
            _gfrd.ReactionRule([self.s1], [self.s1, self.s2]))
        self.assertTrue(True)

        self.m.network_rules.add_reaction_rule(
            _gfrd.ReactionRule([self.s2], [self.s1, self.s2]))
        self.assertTrue(True)

        self.assertRaises(_gfrd.AlreadyExists,
                lambda: self.m.network_rules.add_reaction_rule(
                    _gfrd.ReactionRule([self.s1], [self.s1, self.s2])))

        self.assertRaises(_gfrd.AlreadyExists,
                lambda: self.m.network_rules.add_reaction_rule(
                    _gfrd.ReactionRule([self.s2], [self.s1, self.s2])))

    def test_remove_reaction_rule_1(self):
        # Start with None.
        assert self.m.network_rules.query_reaction_rule(self.s1) == None

        # Add 1.
        rr = _gfrd.ReactionRule([self.s1], [self.s1, self.s2])
        rr['k'] = '0.1'
        self.m.network_rules.add_reaction_rule(rr)

        rules = set(self.m.network_rules.query_reaction_rule(self.s1))
        self.assertEqual(1, len(rules))

        # Remove 1 to get 0.
        self.m.network_rules.remove_reaction_rule(rr)
        gen = self.m.network_rules.query_reaction_rule(self.s1)
        assert len(set(gen)) == 0

    def test_query_reaction_rule(self):
        r1 = _gfrd.ReactionRule([self.s1], [self.s1, self.s2])
        self.m.network_rules.add_reaction_rule(r1)
        a = self.m.network_rules.query_reaction_rule(self.s1)
        self.assertTrue(iter(a) != None)
        a = list(a)
        self.assertEqual(1, len(a))
        self.assertTrue(r1 in a)

        r2 = _gfrd.ReactionRule([self.s1], [self.s1])
        self.m.network_rules.add_reaction_rule(r2)
        a = self.m.network_rules.query_reaction_rule(self.s1)
        self.assertTrue(iter(a) != None)
        a = list(a)
        self.assertEqual(2, len(a))
        self.assertTrue(r1 in a)
        self.assertTrue(r2 in a)

if __name__ == "__main__":
    unittest.main()


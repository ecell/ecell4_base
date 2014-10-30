#!/usr/bin/env python
import _gfrd
import unittest

class NetworkRulesWrapperTestCase(unittest.TestCase):
    def setUp(self):
        self.m = _gfrd.Model()
        self.s1 = _gfrd.SpeciesType()
        self.s1['radius'] = '0.1'
        self.s1['D'] = '0.2'
        self.m.add_species_type(self.s1)
        self.s2 = _gfrd.SpeciesType()
        self.s2['radius'] = '0.3'
        self.s2['D'] = '0.4'
        self.m.add_species_type(self.s2)
        self.nr = _gfrd.NetworkRulesWrapper(self.m.network_rules)

    def tearDown(self):
        pass

    def test(self):
        rr = _gfrd.ReactionRule([self.s1], [self.s1, self.s2])
        rr['k'] = '0.1'
        self.m.network_rules.add_reaction_rule(rr)

        rr = _gfrd.ReactionRule([self.s1, self.s2], [self.s1])
        rr['k'] = '0.2'
        self.m.network_rules.add_reaction_rule(rr)

        rules = set(self.nr.query_reaction_rule(self.s1))
        self.assertEqual(1, len(rules))
        rules = set(self.nr.query_reaction_rule(self.s1, self.s2))
        self.assertEqual(1, len(rules))

if __name__ == "__main__":
    unittest.main()


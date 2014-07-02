from ecell4.core import *

import unittest


class NetowrkModelTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_constructor(self):
        model = NetworkModel()

    def test_add_species_attribute(self):
        model = NetworkModel()
        sp1, sp2 = Species("A"), Species("B")

        self.assertFalse(model.has_species_attribute(sp1))
        self.assertFalse(model.has_species_attribute(sp2))

        model.add_species_attribute(sp1)
        self.assertTrue(model.has_species_attribute(sp1))
        self.assertFalse(model.has_species_attribute(sp2))

        model.remove_species_attribute(sp1)
        self.assertFalse(model.has_species_attribute(sp1))
        self.assertFalse(model.has_species_attribute(sp2))

    def test_query_reaction_rule(self):
        model = NetworkModel()

        sp1, sp2, sp3 = Species("A"), Species("B"), Species("C")
        rr1 = create_degradation_reaction_rule(sp1, 1)
        rr2 = create_unimolecular_reaction_rule(sp1, sp2, 1)
        rr3 = create_binding_reaction_rule(sp1, sp2, sp3, 1)
        rr4 = create_unbinding_reaction_rule(sp3, sp1, sp2, 1)
        model.add_reaction_rule(rr1)
        model.add_reaction_rule(rr2)
        model.add_reaction_rule(rr3)
        model.add_reaction_rule(rr4)
        rules1 = model.query_reaction_rules(sp1)
        rules2 = model.query_reaction_rules(sp2)
        rules3 = model.query_reaction_rules(sp3)
        rules4 = model.query_reaction_rules(sp1, sp2)

        self.assertEqual(len(rules1), 2)
        self.assertEqual(len(rules2), 0)
        self.assertEqual(len(rules3), 1)
        self.assertEqual(len(rules3[0].products()), 2)
        self.assertEqual(len(rules4), 1)
        self.assertEqual(len(rules4[0].products()), 1)
        self.assertEqual(rules4[0].products()[0].serial(), "C")


if __name__ == '__main__':
    unittest.main()

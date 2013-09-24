import unittest

import species
from decorator2 import create_species, create_reaction_rule


class ReactionRuleTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test_creation(self):
        sp1 = create_species("A")
        sp2 = create_species("B")
        sp3 = create_species("C")
        rr1 = species.ReactionRule([sp1, sp2], [sp3])

        self.assertEqual(str(rr1), "A+B>C")
        self.assertEqual(str(rr1), str(create_reaction_rule("A+B>C")))
        # self.assertEqual(rr1, create_reaction_rule("A+B>C"))

    def test_properties1(self):
        rr1 = create_reaction_rule("A+B>C")

        self.assertEqual(rr1.num_reactants(), 2)
        self.assertEqual(rr1.num_products(), 1)
        self.assertEqual(rr1.options(), [])
        self.assertEqual(rr1.reactants()[0], create_species("A"))
        self.assertEqual(rr1.reactants()[1], create_species("B"))
        self.assertEqual(rr1.products()[0], create_species("C"))

    def test_properties2(self):
        rr1 = create_reaction_rule("A+B>C")
        self.assertFalse(rr1.is_degradation() or rr1.is_synthesis())

        rr2 = create_reaction_rule("A>~A")
        self.assertTrue(rr2.is_degradation())
        self.assertFalse(rr2.is_synthesis())

        rr3 = create_reaction_rule("~A>A")
        self.assertFalse(rr3.is_degradation())
        self.assertTrue(rr3.is_synthesis())

    def test_matches1(self):
        rr1 = create_reaction_rule("A+B>C")

        self.assertEqual(len(rr1.match(create_species("B"), create_species("A"))), 0)
        self.assertEqual(len(rr1.match(create_species("A"))), 0)
        self.assertEqual(len(rr1.match(
            create_species("A"), create_species("B"), create_species("C"))), 0)

        self.assertEqual(len(rr1.match(create_species("A"), create_species("B"))), 1)
        self.assertEqual(len(rr1.match(create_species("A"), create_species("B"))[0]), 1)
        self.assertEqual(
            rr1.match(create_species("A"), create_species("B"))[0][0],
            create_species("C"))


if __name__ == '__main__':
    unittest.main()

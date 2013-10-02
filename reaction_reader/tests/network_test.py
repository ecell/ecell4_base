import unittest

import network
from decorator2 import create_species, create_reaction_rule


def dump_reactions(reactions):
    for i, reaction in enumerate(reactions):
        print "[%d] %s > %s | %s" % (
            i + 1, "+".join([str(elem) for elem in reaction[0]]),
            "+".join([str(elem) for elem in reaction[1]]),
            reaction[2])

class NetworkTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test_check_stoichiometry(self):
        sp1 = create_species(
            "A(l,r^1).A(l^1,r^2).B(l^2,r^3).A(l^3,r^4)"
            + ".B(l^4,r^5).A(l^5,r^6).B(l^6,r)")
        self.assertTrue(
            network.check_stoichiometry(sp1, dict(A=4, B=3)))
        self.assertTrue(
            network.check_stoichiometry(sp1, dict(C=0)))
        self.assertFalse(
            network.check_stoichiometry(sp1, dict(A=3, B=3)))
        self.assertFalse(
            network.check_stoichiometry(sp1, dict(A=4, B=2)))

    def test_generate_recurse1(self):
        rr1 = create_reaction_rule("X(c=u)>X(c=p)")
        rr2 = create_reaction_rule("X(a)+X(a)>X(a^1).X(a^1)")
        sp1 = create_species("X(a,b,c=u)")

        newseeds, seeds, newreactions = network.generate_recurse2(
            [sp1], [rr1], [], {})
        self.assertEqual(len(seeds), 1)
        self.assertEqual(seeds[0], sp1)
        self.assertEqual(len(newseeds), 1)
        self.assertEqual(newseeds[0], create_species("X(a,b,c=p)"))
        self.assertEqual(len(newreactions), 1)
        # dump_reactions(newreactions)

        newseeds, seeds, newreactions = network.generate_recurse2(
            [sp1], [rr1, rr2], [], {})
        self.assertEqual(len(seeds), 1)
        self.assertEqual(seeds[0], sp1)
        self.assertEqual(len(newseeds), 2)
        self.assertIn(create_species("X(a,b,c=p)"), newseeds)
        self.assertIn(create_species("X(a^1,b,c=u).X(a^1,b,c=u)"), newseeds)
        self.assertEqual(len(newreactions), 2)
        # dump_reactions(newreactions)

        newseeds, seeds, newreactions = network.generate_recurse2(
            newseeds, [rr1, rr2], seeds, {})
        self.assertEqual(len(seeds), 3)
        self.assertEqual(len(newseeds), 2)
        self.assertIn(create_species("X(a^1,b,c=p).X(a^1,b,c=u)"), newseeds)
        self.assertIn(create_species("X(a^1,b,c=p).X(a^1,b,c=p)"), newseeds)
        self.assertEqual(len(newreactions), 5)
        # dump_reactions(newreactions)

    def test_generate_recurse2(self):
        sp1, sp2 = create_species("A(l,r)"), create_species("B(l)")

        rr1 = create_reaction_rule("A(r)+B(l)>A(r^1).B(l^1)")
        newseeds, seeds, newreactions = network.generate_recurse2(
            [sp1, sp2], [rr1], [], {})
        self.assertEqual(len(newseeds), 1)
        self.assertEqual(newseeds[0], create_species("A(l,r^1).B(l^1)"))
        self.assertEqual(len(newreactions), 1)

        rr2 = create_reaction_rule("A(r)+A(l,r)+B(l)>A(r^1).A(l^1,r^2).B(l^2)")
        newseeds, seeds, newreactions = network.generate_recurse2(
            [sp1, sp2], [rr2], [], {})
        self.assertEqual(len(newseeds), 1)
        self.assertEqual(
            newseeds[0], create_species("A(l,r^1).A(l^1,r^2).B(l^2)"))
        self.assertEqual(len(newreactions), 1)
        # dump_reactions(newreactions)


if __name__ == '__main__':
    unittest.main()

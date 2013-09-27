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

    def test_matches2(self):
        sp1 = create_species("X(a=b1,c=d^1).Y(e^1)")

        retval = create_reaction_rule("X(a=b1)>X(a=b2)").match(sp1)
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 1)
        self.assertEqual(retval[0][0], create_species("X(a=b2,c=d^1).Y(e^1)"))

        retval = create_reaction_rule("X(a=b1)>X(a=b2)+Y(a=b)").match(sp1)
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 2)
        self.assertTrue(create_species("X(a=b2,c=d^1).Y(e^1)") in retval[0])
        self.assertTrue(create_species("Y(a=b)") in retval[0])

        retval = create_reaction_rule("X(a=b1)>~X(a=b1)").match(sp1)
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 0)

        retval = create_reaction_rule("X(c^1).Y(e^1)>X(c)+Y(e)").match(sp1)
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 2)
        self.assertTrue(create_species("Y(e)") in retval[0])
        self.assertTrue(create_species("X(a=b1,c=d)") in retval[0])

        retval = create_reaction_rule("X(c^1).Y(e^1)>~X(c)+Y(e)").match(sp1)
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 1)
        self.assertEqual(retval[0][0], create_species("Y(e)"))

        retval = create_reaction_rule("X(a)+~Z(a)>X(a^1).Z(a^1)").match(sp1)
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 1)
        self.assertEqual(
            retval[0][0], create_species("X(a=b1^2,c=d^1).Y(e^1).Z(a^2)"))

        retval = create_reaction_rule("X(a=b1)>X(a=b2)").match(
            create_species("X(a=b1,c=d1^1).X(a=b1,c=d2^1)"))
        self.assertEqual(len(retval), 2)
        self.assertEqual(len(retval[0]), 1)
        self.assertEqual(len(retval[1]), 1)
        products = [elem[0] for elem in retval]
        self.assertTrue(
            create_species("X(a=b2,c=d1^1).X(a=b1,c=d2^1)") in products)
        self.assertTrue(
            create_species("X(a=b1,c=d1^1).X(a=b2,c=d2^1)") in products)

    def test_matches3(self):
        sp1 = create_species("X(l,r)")
        sp2 = create_species("Y(l,r^1).Z(l^1,r)")

        retval = create_reaction_rule("X(r)+Y(l)>X(r^1).Y(l^1)").match(sp1, sp2)
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 1)
        self.assertEqual(
            retval[0][0], create_species("X(l,r^1).Y(l^1,r^2).Z(l^2,r)"))

        retval = create_reaction_rule("X(r)+Y(r)>X(r^1).Y(r^1)").match(sp1, sp2)
        self.assertEqual(len(retval), 0)

        retval = create_reaction_rule("X(r)+Y(l,r^_)>X(r^1).Y(l^1,r^_)").match(sp1, sp2)
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 1)
        self.assertEqual(
            retval[0][0], create_species("X(l,r^1).Y(l^1,r^2).Z(l^2,r)"))

        retval = create_reaction_rule("X(r)+Y(l,r)>X(r^1).Y(l^1,r)").match(sp1, sp2)
        self.assertEqual(len(retval), 0)

        retval = create_reaction_rule("X(r)+Y(l)>X(r^1).Y(l^1)").match(sp2, sp1)
        self.assertEqual(len(retval), 0)

    def test_matches4(self):
        sp1 = create_species("X(a=b)")
        sp2 = create_species("Y(a=b,c=d)")

        retval = create_reaction_rule("_(a=b)>_(a=c)").match(sp1)
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 1)
        self.assertEqual(retval[0][0], create_species("X(a=c)"))

        retval = create_reaction_rule("_(a=b)>_(a=c)").match(
            create_species("X(c^1,a=b).Y(c^1,d^2,a=b).X(c^2,a=b)"))
        self.assertEqual(len(retval), 3)

        retval = create_reaction_rule("_1(a=b)>_1(a=c)").match(sp1)
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 1)
        self.assertEqual(retval[0][0], create_species("X(a=c)"))

        retval = create_reaction_rule("_1(a)+_1(a)>_1(a^1)._1(a^1)").match(
            sp1, sp1)
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 1)
        self.assertEqual(retval[0][0], create_species("X(a=b^1).X(a=b^1)"))

        retval = create_reaction_rule("_1(a)+_1(a)>_1(a^1)._1(a^1)").match(
            sp1, sp2)
        self.assertEqual(len(retval), 0)

        retval = create_reaction_rule("_1(a)+_2(a)>_1(a^1)._2(a^1)").match(
            sp1, sp1)
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 1)
        self.assertEqual(retval[0][0], create_species("X(a=b^1).X(a=b^1)"))

        retval = create_reaction_rule("_1(a)+_2(a)>_1(a^1)._2(a^1)").match(
            sp1, sp2)
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 1)
        self.assertEqual(retval[0][0], create_species("X(a=b^1).Y(a=b^1,c=d)"))

    def test_matches5(self):
        sp1 = create_species("X(a=c13,b=c12)")
        sp2 = create_species("Y(a=c12,b=c13)")

        retval = create_reaction_rule("X(a=_1,b=_2)>Y(a=_2,b=_1)").match(sp1)
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 1)
        self.assertEqual(retval[0][0], sp2)

        retval = create_reaction_rule(
            "X(a=_1)+Y(b=_1)>X(a=_1^1).Y(b=_1^1)").match(sp1, sp2)
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 1)
        self.assertEqual(
            retval[0][0], create_species("X(a=c13^1,b=c12).Y(a=c12,b=c13^1)"))

        retval = create_reaction_rule(
            "X(a=_1)+Y(a=_1)>X(a=_1^1).Y(a=_1^1)").match(sp1, sp2)
        self.assertEqual(len(retval), 0)

        retval = create_reaction_rule(
            "X(a=_1)+Y(a=_2)>X(a=_2^1).Y(a=_1^1)").match(sp1, sp2)
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 1)
        self.assertEqual(
            retval[0][0], create_species("X(a=c12^1,b=c12).Y(a=c13^1,b=c13)"))

    def test_matches6(self):
        sp1 = create_species("X(a,b,c)")

        retval = create_reaction_rule("X>Y").match(sp1)
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 1)
        self.assertEqual(retval[0][0], create_species("Y"))

        retval = create_reaction_rule("X(c)>Y").match(sp1)
        self.assertEqual(len(retval), 1)

        retval = create_reaction_rule("X(~c)>Y").match(sp1)
        self.assertEqual(len(retval), 0)

        retval = create_reaction_rule("X(a,~d)>Y").match(sp1)
        self.assertEqual(len(retval), 1)


if __name__ == '__main__':
    unittest.main()

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
        sp1 = create_species("A")

        self.assertEqual(len(rr1.generate([create_species("B"), sp1])), 0)
        self.assertEqual(len(rr1.generate([sp1])), 0)
        self.assertEqual(len(rr1.generate(
            [sp1, create_species("B"), create_species("C")])), 0)

        retval = rr1.generate([sp1, create_species("B")])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 1)
        self.assertEqual(
            retval[0][1][0], create_species("C"))

        rr2 = create_reaction_rule("A>A+A")
        retval = rr2.generate([sp1])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 2)
        self.assertEqual(retval[0][1][0], retval[0][1][1])
        self.assertEqual(retval[0][1][0], sp1)

    def test_matches2(self):
        sp1 = create_species("X(a=b1,c=d^1).Y(e^1)")

        retval = create_reaction_rule("X(a=b1)>X(a=b2)").generate([sp1])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 1)
        self.assertEqual(retval[0][1][0], create_species("X(a=b2,c=d^1).Y(e^1)"))

        retval = create_reaction_rule("X(a=b1)>X(a=b2)+Y(a=b)").generate([sp1])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 2)
        self.assertTrue(create_species("X(a=b2,c=d^1).Y(e^1)") in retval[0][1])
        self.assertTrue(create_species("Y(a=b)") in retval[0][1])

        retval = create_reaction_rule("X(c^1).Y(e^1)>X(c)+Y(e)").generate([sp1])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 2)
        self.assertTrue(create_species("Y(e)") in retval[0][1])
        self.assertTrue(create_species("X(a=b1,c=d)") in retval[0][1])

        retval = create_reaction_rule("X(c^1).Y(e^1)>~X(c)+Y(e)").generate(
            [sp1])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 1)
        self.assertEqual(retval[0][1][0], create_species("Y(e)"))

        retval = create_reaction_rule("X(a)+~Z(a)>X(a^1).Z(a^1)").generate([sp1])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 1)
        self.assertEqual(
            retval[0][1][0], create_species("X(a=b1^2,c=d^1).Y(e^1).Z(a^2)"))

        retval = create_reaction_rule("X(a=b1)>X(a=b2)").generate(
            [create_species("X(a=b1,c=d1^1).X(a=b1,c=d2^1)")])
        self.assertEqual(len(retval), 2)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[1]), 3)
        self.assertEqual(len(retval[0][1]), 1)
        self.assertEqual(len(retval[1][1]), 1)
        products = [elem[1][0] for elem in retval]
        self.assertTrue(
            create_species("X(a=b2,c=d1^1).X(a=b1,c=d2^1)") in products)
        self.assertTrue(
            create_species("X(a=b1,c=d1^1).X(a=b2,c=d2^1)") in products)


        retval = create_reaction_rule("X(c^_,a=b1)>X(c^_,a=b2)").generate([sp1])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 1)
        self.assertEqual(retval[0][1][0], create_species("X(a=b2,c=d^1).Y(e^1)"))

    def test_matches3(self):
        sp1 = create_species("X(l,r)")
        sp2 = create_species("Y(l,r^1).Z(l^1,r)")

        retval = create_reaction_rule("X(r)+Y(l)>X(r^1).Y(l^1)").generate(
            [sp1, sp2])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 1)
        self.assertEqual(
            retval[0][1][0], create_species("X(l,r^1).Y(l^1,r^2).Z(l^2,r)"))

        retval = create_reaction_rule("X(r)+Y(r)>X(r^1).Y(r^1)").generate(
            [sp1, sp2])
        self.assertEqual(len(retval), 0)

        retval = create_reaction_rule("X(r)+Y(l,r^_)>X(r^1).Y(l^1,r^_)").generate(
            [sp1, sp2])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 1)
        self.assertEqual(
            retval[0][1][0], create_species("X(l,r^1).Y(l^1,r^2).Z(l^2,r)"))

        retval = create_reaction_rule("X(r)+Y(l,r)>X(r^1).Y(l^1,r)").generate(
            [sp1, sp2])
        self.assertEqual(len(retval), 0)

        retval = create_reaction_rule("X(r)+Y(l)>X(r^1).Y(l^1)").generate(
            [sp2, sp1])
        self.assertEqual(len(retval), 0)

    def test_matches4(self):
        sp1 = create_species("X(a=b)")
        sp2 = create_species("Y(a=b,c=d)")

        retval = create_reaction_rule("_(a=b)>_(a=c)").generate([sp1])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 1)
        self.assertEqual(retval[0][1][0], create_species("X(a=c)"))

        retval = create_reaction_rule("_(a=b)>_(a=c)").generate(
            [create_species("X(c^1,a=b).Y(c^1,d^2,a=b).X(c^2,a=b)")])
        self.assertEqual(len(retval), 3)

        retval = create_reaction_rule("_1(a=b)>_1(a=c)").generate([sp1])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 1)
        self.assertEqual(retval[0][1][0], create_species("X(a=c)"))

        retval = create_reaction_rule("_1(a)+_1(a)>_1(a^1)._1(a^1)").generate(
            [sp1, sp1])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 1)
        self.assertEqual(retval[0][1][0], create_species("X(a=b^1).X(a=b^1)"))

        retval = create_reaction_rule("_1(a)+_1(a)>_1(a^1)._1(a^1)").generate(
            [sp1, sp2])
        self.assertEqual(len(retval), 0)

        retval = create_reaction_rule("_1(a)+_2(a)>_1(a^1)._2(a^1)").generate(
            [sp1, sp1])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 1)
        self.assertEqual(retval[0][1][0], create_species("X(a=b^1).X(a=b^1)"))

        retval = create_reaction_rule("_1(a)+_2(a)>_1(a^1)._2(a^1)").generate(
            [sp1, sp2])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 1)
        self.assertEqual(
            retval[0][1][0], create_species("X(a=b^1).Y(a=b^1,c=d)"))

    def test_matches5(self):
        sp1 = create_species("X(a=c13,b=c12)")
        sp2 = create_species("Y(a=c12,b=c13)")

        retval = create_reaction_rule("X(a=_1,b=_2)>Y(a=_2,b=_1)").generate(
            [sp1])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 1)
        self.assertEqual(retval[0][1][0], sp2)

        retval = create_reaction_rule(
            "X(a=_1)+Y(b=_1)>X(a=_1^1).Y(b=_1^1)").generate([sp1, sp2])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 1)
        self.assertEqual(
            retval[0][1][0],
            create_species("X(a=c13^1,b=c12).Y(a=c12,b=c13^1)"))

        retval = create_reaction_rule(
            "X(a=_1)+Y(a=_1)>X(a=_1^1).Y(a=_1^1)").generate([sp1, sp2])
        self.assertEqual(len(retval), 0)

        retval = create_reaction_rule(
            "X(a=_1)+Y(a=_2)>X(a=_2^1).Y(a=_1^1)").generate([sp1, sp2])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 1)
        self.assertEqual(
            retval[0][1][0],
            create_species("X(a=c12^1,b=c12).Y(a=c13^1,b=c13)"))

    def test_matches6(self):
        sp1 = create_species("X(a,b,c)")

        retval = create_reaction_rule("X>Y").generate([sp1])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 1)
        self.assertEqual(retval[0][1][0], create_species("Y"))

        retval = create_reaction_rule("X(c)>Y").generate([sp1])
        self.assertEqual(len(retval), 1)

        retval = create_reaction_rule("X(~c)>Y").generate([sp1])
        self.assertEqual(len(retval), 0)

        retval = create_reaction_rule("X(a,~d)>Y").generate([sp1])
        self.assertEqual(len(retval), 1)

    def test_matches7(self):
        rr1 = create_reaction_rule("_(ps=u)>_(ps=p)")
        self.assertEqual(
            len(rr1.generate([create_species("A(ps1=u,ps2=u)")])), 0)

        sp1 = create_species("A(ps1=u,ps2=u,ps=[ps1,ps2])")
        retval = rr1.generate([sp1])
        self.assertEqual(len(retval), 2)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[1]), 3)
        self.assertEqual(len(retval[0][1]), 1)
        self.assertEqual(len(retval[1][1]), 1)
        self.assertTrue(
            create_species("A(ps1=p,ps2=u,ps=[ps1,ps2])")
            in (retval[0][1][0], retval[1][1][0]))
        self.assertTrue(
            create_species("A(ps1=u,ps2=p,ps=[ps1,ps2])")
            in (retval[0][1][0], retval[1][1][0]))

        rr2 = create_reaction_rule("_(ps)+_(ps)>_(ps^1)._(ps^1)")
        self.assertEqual(len(rr2.generate([sp1, sp1])), 4)

    def test_matches8(self):
        rr1 = create_reaction_rule(
            "A(_1=u,_2=u,ps=[_1,_2])>A(_1=p,_2=u,ps=[_1,_2])")
        sp1 = create_species("A(ps1=u,ps2=u,ps=[ps1,ps2])")

        self.assertEqual(
            len(rr1.generate([create_species("A(ps1=u,ps=[ps1])")])), 0)
        self.assertEqual(len(rr1.generate([sp1])), 2)
        self.assertEqual(
            len(rr1.generate(
                [create_species("A(ps1=u,ps2=u,ps3=u,ps=[ps1,ps2,ps3])")])), 6)

        rr2 = create_reaction_rule(
            "A(_1=u,ps=[_1]) + A(_1=u,ps=[_1])"
            " > A(_1=u^1,ps=[_1]).A(_1=u^1,ps=[_1])")
        rr3 = create_reaction_rule(
            "A(_1=u,ps=[_1]) + A(_2=u,ps=[_2])"
            " > A(_1=u^1,ps=[_1]).A(_2=u^1,ps=[_2])")
        self.assertEqual(len(rr2.generate([sp1, sp1])), 2)
        self.assertEqual(len(rr3.generate([sp1, sp1])), 4)

    def test_matches9(self):
        rr1 = create_reaction_rule(
            "A(r^1).B(bs^1)+C(bs)+D(l)>A(r^1).C(bs^1)+B(bs^1).D(l^1)")
        sp1, sp2, sp3, sp4 = (
            create_species("A(l,r)"),
            create_species("B(bs)"),
            create_species("C(bs)"),
            create_species("D(l,r)"))
        sp5 = create_species("A(l^2,r^1).B(bs^1).C(bs^2)")
        sp6 = create_species("D(l,r^1).B(bs^1)")

        retval = rr1.generate([sp1, sp3, sp4])
        self.assertEqual(len(retval), 0)

        retval = rr1.generate([sp5, sp3, sp6])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 2)
        self.assertIn(create_species("A(l^1,r^2).C(bs^1).C(bs^2)"), retval[0][1])
        self.assertIn(create_species("B(bs^1).D(l^1,r^2).B(bs^2)"), retval[0][1])

        rr2 = create_reaction_rule("A+B+C>D+E+F")
        retval = rr2.generate(
            [create_species("A"), create_species("B"), create_species("C")])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 3)
        self.assertIn(create_species("D"), retval[0][1])
        self.assertIn(create_species("E"), retval[0][1])
        self.assertIn(create_species("F"), retval[0][1])

    def test_matches10(self):
        rr1 = create_reaction_rule("R(l) + R(l) > R(l^1).R(l^1)")
        rr2 = create_reaction_rule(
            "R(_1,_2,l=[_1,_2]) + R(l) > R(_1^1,_2,l=[_1,_2]).R(l^1)")
        rr3 = create_reaction_rule("R(l1,l2) + R(l) > R(l1^1,l2).R(l^1)")

        sp1 = create_species("R(l1,l2,l=(l1,l2))")

        for rr in (rr1, rr2, rr3):
            retval = rr.generate([sp1, sp1])
            self.assertEqual(len(retval), 4)
            self.assertTrue(
                all([len(elem) == 3 and len(elem[1]) == 1 for elem in retval]))
            for i in range(4):
                for j in range(i, 4):
                    self.assertEqual(retval[i][1][0], retval[j][1][0])

        rr4 = create_reaction_rule(
            "R(_1,l=[_1]) + R(_1,l=[_1]) > R(_1^1,l=[_1]).R(_1^1,l=[_1])")
        retval = rr4.generate([sp1, sp1])
        self.assertEqual(len(retval), 2)
        self.assertTrue(
            all([len(elem) == 3 and len(elem[1]) == 1 for elem in retval]))
        self.assertEqual(retval[0][1][0], retval[1][1][0])

    def test_exceptions1(self):
        rr1 = create_reaction_rule("A(_1=u)>A(_1=p)")
        sp1 = create_species("A(ps1=u,ps2=u)")
        self.assertRaises(RuntimeError, rr1.generate, [sp1])

        self.assertRaises(
            RuntimeError, create_reaction_rule,
            "A(ps1^_1).A(ps1^_1)>A(ps1)+A(ps1)")

    def test_options1(self):
        rr1 = create_reaction_rule("A(r)+B(l)>A(r^1).B(l^1)|ExcludeReactants(1,B)")
        sp1 = create_species("A(l,r)")
        sp2 = create_species("B(l,r)")
        self.assertEqual(len(rr1.generate([sp1, sp2])), 1)
        self.assertEqual(
            len(rr1.generate([create_species("A(l^1,r).B(l,r^1)"), sp2])), 0)
        self.assertEqual(
            len(rr1.generate([create_species("A(l^1,r).C(l,r^1)"), sp2])), 1)

        rr2 = create_reaction_rule("A(r)+B(l)>A(r^1).B(l^1)|IncludeReactants(1,B)")
        self.assertEqual(len(rr2.generate([sp1, sp2])), 0)
        self.assertEqual(
            len(rr2.generate([create_species("A(l^1,r).B(l,r^1)"), sp2])), 1)
        self.assertEqual(
            len(rr2.generate([create_species("A(l^1,r).C(l,r^1)"), sp2])), 0)

        rr3 = create_reaction_rule("A(r)+B(l)>A(r^1).B(l^1)|ExcludeProducts(1,B)")
        self.assertEqual(len(rr3.generate([sp1, sp2])), 0)
        self.assertEqual(
            len(rr3.generate([create_species("A(l^1,r).B(l,r^1)"), sp2])), 0)
        self.assertEqual(
            len(rr3.generate([create_species("A(l^1,r).C(l,r^1)"), sp2])), 0)

        rr4 = create_reaction_rule("A(r)+B(l)>A(r^1).B(l^1)|IncludeProducts(1,B)")
        self.assertEqual(len(rr4.generate([sp1, sp2])), 1)
        self.assertEqual(
            len(rr4.generate([create_species("A(l^1,r).B(l,r^1)"), sp2])), 1)
        self.assertEqual(
            len(rr4.generate([create_species("A(l^1,r).C(l,r^1)"), sp2])), 1)

        rr5 = create_reaction_rule("_1(r) + _2(l) > _1(r^1)._2(l^1)|CaseIf(1.0,_1=A)")
        retval = rr5.generate([sp1, sp2])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(retval[0][2], 1.0)

        retval = rr5.generate([sp2, sp1])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(retval[0][2], 0.0) # default

    def test_failures(self):
        return # exit with no test

        rr1 = create_reaction_rule("A>B")
        self.assertEqual(len(rr1.generate([create_species("A(a^1).B(b^1)")])), 0)

        sp1 = create_species("X(a=b1,c=d^1).Y(e^1)")
        retval = create_reaction_rule("X(a=b1)>~X(a=b1)").generate([sp1])
        self.assertEqual(len(retval), 1)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 0)


if __name__ == '__main__':
    unittest.main()

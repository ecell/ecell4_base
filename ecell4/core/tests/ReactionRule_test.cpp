#define BOOST_TEST_MODULE "ReactionRule_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/floating_point_comparison.hpp>

#include <ecell4/core/Species.hpp>
#include <ecell4/core/Context.hpp>
#include <ecell4/core/ReactionRule.hpp>

using namespace ecell4;


BOOST_AUTO_TEST_CASE(ReactionRule_test_constructor)
{
    ReactionRule rr;
}

BOOST_AUTO_TEST_CASE(ReactionRule_test_k)
{
    const Real epsrel(1e-6);
    const Real k(1.5);
    ReactionRule rr;
    rr.set_k(0);
    rr.set_k(k);
    BOOST_CHECK_CLOSE(rr.k(), k, k * epsrel);
    BOOST_CHECK_THROW(rr.set_k(-k), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(ReactionRule_test_reactants)
{
    ReactionRule rr;
    Species sp1("A"), sp2("B");
    rr.add_reactant(sp1);
    rr.add_reactant(sp2);
    rr.add_reactant(sp1);

    const ReactionRule::reactant_container_type& reactants(rr.reactants());
    BOOST_CHECK_EQUAL(reactants.size(), 3);
    BOOST_CHECK_EQUAL(std::count(reactants.begin(), reactants.end(), sp1), 2);
    BOOST_CHECK_EQUAL(std::count(reactants.begin(), reactants.end(), sp2), 1);
}

BOOST_AUTO_TEST_CASE(ReactionRule_test_products)
{
    ReactionRule rr;
    Species sp1("A"), sp2("B");
    rr.add_product(sp1);
    rr.add_product(sp2);
    rr.add_product(sp1);

    const ReactionRule::product_container_type& products(rr.products());
    BOOST_CHECK_EQUAL(products.size(), 3);
    BOOST_CHECK_EQUAL(std::count(products.begin(), products.end(), sp1), 2);
    BOOST_CHECK_EQUAL(std::count(products.begin(), products.end(), sp2), 1);
}

BOOST_AUTO_TEST_CASE(ReactionRule_test_compare)
{
    ReactionRule rr1, rr2, rr3, rr4, rr5;
    Species sp1("A"), sp2("B"), sp3("C"), sp4("D");

    rr1.add_reactant(sp1);
    rr1.add_reactant(sp2);
    rr1.add_product(sp3);
    rr1.set_k(1.5);

    rr2.add_reactant(sp1);
    rr2.add_reactant(sp3);
    rr2.add_product(sp3);
    rr2.set_k(1.5);

    rr3.add_reactant(sp1);
    rr3.add_reactant(sp2);
    rr3.add_product(sp4);
    rr3.set_k(1.5);

    rr4.add_reactant(sp2);
    rr4.add_reactant(sp1);
    rr4.add_product(sp3);
    rr4.set_k(5.0);

    rr5.add_reactant(sp1);
    rr5.add_reactant(sp2);
    rr5.add_product(sp3);
    rr5.set_k(5.0);

    BOOST_CHECK(rr1 == rr1);
    BOOST_CHECK(rr1 != rr2);
    BOOST_CHECK(rr1 != rr3);
    BOOST_CHECK(rr1 != rr4);
    BOOST_CHECK(rr1 == rr5);
}

BOOST_AUTO_TEST_CASE(ReactionRule_test_generate1)
{
    ReactionRule rr1;
    Species sp1("A");
    rr1.add_product(sp1);
    rr1.set_k(1.0);

    std::vector<Species> reactants;  // XXX: empty
    std::vector<ReactionRule> retval(rr1.generate(reactants));

    BOOST_CHECK_EQUAL(retval.size(), 1);
    BOOST_CHECK(retval[0] == rr1);
}

BOOST_AUTO_TEST_CASE(ReactionRule_test_generate2)
{
    ReactionRule rr1;
    rr1.add_reactant(Species("_1(b)"));
    rr1.add_reactant(Species("_1(b)"));
    rr1.add_product(Species("_1(b^1)._1(b^1)"));

    ReactionRule::reactant_container_type reactants1;
    reactants1.push_back(Species("A(a^1,b).B(a^1,b)"));
    reactants1.push_back(Species("B(a,b)"));

    BOOST_CHECK_EQUAL(rr1.count(reactants1), 1);

    std::vector<ReactionRule> retval;
    retval = rr1.generate(reactants1);
    BOOST_CHECK_EQUAL(retval.size(), 1);
    BOOST_CHECK_EQUAL(retval[0].products().size(), 1);
    BOOST_CHECK_EQUAL(retval[0].products()[0].units().size(), 3);

    ReactionRule rr2;
    rr2.add_reactant(Species("A(b)"));
    rr2.add_reactant(Species("B(b)"));
    rr2.add_product(Species("A(b^1).B(b^1)"));
    ReactionRule::reactant_container_type reactants2;
    reactants2.push_back(Species("A(a^1,b).A(a^1,b)"));
    reactants2.push_back(Species("B(a^1,b).B(a^1,b^2).B(a^2,b)"));

    BOOST_CHECK_EQUAL(rr2.count(reactants2), 4);

    ReactionRule rr3;
    rr3.add_reactant(Species("A"));
    rr3.add_product(Species("B"));
    ReactionRule::reactant_container_type reactants3;
    reactants3.push_back(Species("A"));
    retval = rr3.generate(reactants3);
    BOOST_CHECK_EQUAL(retval.size(), 1);
    BOOST_CHECK_EQUAL(retval[0].products().size(), 1);
    BOOST_CHECK_EQUAL(retval[0].products()[0].serial(), "B");

    ReactionRule rr4;
    rr4.add_reactant(Species("A(b^1).B(b^1)"));
    rr4.add_product(Species("A(b)"));
    rr4.add_product(Species("B(b)"));
    ReactionRule::reactant_container_type reactants4;
    reactants4.push_back(
        Species("A(a^1,b^5).A(a^1,b^4).B(a^2,b).B(a^2,b^3).B(a^3,b^4).B(a,b^5)"));
    retval = rr4.generate(reactants4);
    BOOST_CHECK_EQUAL(retval.size(), 2);
    BOOST_CHECK_EQUAL(retval[0].products().size(), 2);
    BOOST_CHECK_EQUAL(retval[1].products().size(), 2);
    BOOST_CHECK_EQUAL(retval[0].products()[0].units().size() + retval[0].products()[1].units().size(), 6);
    BOOST_CHECK_EQUAL(retval[1].products()[0].units().size() + retval[1].products()[1].units().size(), 6);
}

BOOST_AUTO_TEST_CASE(ReactionRule_test_recursive_generation1)
{
    ReactionRule rr1;
    rr1.add_reactant(Species("X(r^1).X(l^1)"));
    rr1.add_product(Species("X(r)"));
    rr1.add_product(Species("X(l)"));
    rr1.set_k(1.0);
    const Species sp1("X(l,r^1).X(l^1,r^2).X(l^2,r^3).X(l^3,r^4).X(l^4,r)");

    ReactionRuleExpressionMatcher rrexp(rr1);
    BOOST_CHECK(rrexp.match(sp1));

    unsigned int i(0);
    do {
        ++i;
        std::vector<Species> products(rrexp.generate());
        // const ReactionRule tmp(rrexp.reactants(), products, rr1.k());
        // std::cerr << "GEN: " << tmp.as_string() << std::endl;
    } while (rrexp.next());

    BOOST_CHECK_EQUAL(i, 4);

    std::vector<ReactionRule> retval(rr1.generate(rrexp.reactants()));
    BOOST_CHECK_EQUAL(retval.size(), 4);
}

BOOST_AUTO_TEST_CASE(ReactionRule_test_generate3)
{
    ReactionRule rr1;
    rr1.add_reactant(Species("A(b=u^1).A(b=u^1)"));
    rr1.add_product(Species("A(b=u)"));
    rr1.add_product(Species("A(b=u)"));
    rr1.set_k(1.0);

    ReactionRule rr2;
    rr2.add_reactant(Species("A(b=u^1).A(b=u^1)"));
    rr2.add_product(Species("A(b=u)"));
    rr2.add_product(Species("A(b=p)"));
    rr2.set_k(1.0);

    ReactionRule::reactant_container_type reactants1(1, Species("A(b=u^1).A(b=u^1)"));

    std::vector<ReactionRule> retval1 = rr1.generate(reactants1);
    BOOST_CHECK_EQUAL(retval1.size(), 1);
    BOOST_CHECK_EQUAL(retval1[0].k(), 1.0);
    BOOST_CHECK_EQUAL(retval1[0].reactants().size(), 1);
    BOOST_CHECK_EQUAL(retval1[0].reactants()[0], reactants1[0]);
    BOOST_CHECK_EQUAL(retval1[0].products().size(), 2);
    BOOST_CHECK_EQUAL(retval1[0].products()[0], Species("A(b=u)"));
    BOOST_CHECK_EQUAL(retval1[0].products()[1], Species("A(b=u)"));

    std::vector<ReactionRule> retval2 = rr2.generate(reactants1);
    BOOST_CHECK_EQUAL(retval2.size(), 2);
    BOOST_CHECK_EQUAL(retval2[0].k(), 1.0);
    BOOST_CHECK_EQUAL(retval2[1].k(), 1.0);
    BOOST_CHECK_EQUAL(retval2[0].reactants().size(), 1);
    BOOST_CHECK_EQUAL(retval2[1].reactants().size(), 1);
    BOOST_CHECK_EQUAL(retval2[0].reactants()[0], reactants1[0]);
    BOOST_CHECK_EQUAL(retval2[1].reactants()[0], reactants1[0]);
    BOOST_CHECK_EQUAL(retval2[0].products().size(), 2);
    BOOST_CHECK_EQUAL(retval2[1].products().size(), 2);
    BOOST_CHECK(
        (retval2[0].products()[0] == Species("A(b=u)")
         && retval2[0].products()[1] == Species("A(b=p)"))
        || (retval2[0].products()[0] == Species("A(b=p)")
         && retval2[0].products()[1] == Species("A(b=u)")));
    BOOST_CHECK(
        (retval2[1].products()[0] == Species("A(b=u)")
         && retval2[1].products()[1] == Species("A(b=p)"))
        || (retval2[1].products()[0] == Species("A(b=p)")
         && retval2[1].products()[1] == Species("A(b=u)")));

    ReactionRule rr3;
    rr3.add_reactant(Species("A(b^1,c^2).A(b^1,c^3).B(l^2,r^4).B(l^3,r^4)"));
    rr3.add_product(Species("A(b^1,c).A(b^1,c)"));
    rr3.add_product(Species("B(l,r^1).B(l,r^1)"));
    rr3.set_k(1.0);

    ReactionRule::reactant_container_type reactants2(1, Species("A(b^1,c^2).A(b^1,c^3).B(l^2,r^4).B(l^3,r^4)"));
    std::vector<ReactionRule> retval3 = rr3.generate(reactants2);

    BOOST_CHECK_EQUAL(retval3.size(), 1);
    BOOST_CHECK_EQUAL(retval3[0].k(), 1.0);
    BOOST_CHECK_EQUAL(retval3[0].reactants().size(), 1);
    BOOST_CHECK_EQUAL(retval3[0].reactants()[0], reactants2[0]);
    BOOST_CHECK_EQUAL(retval3[0].products().size(), 2);
    BOOST_CHECK_EQUAL(retval3[0].products()[0], Species("A(b^1,c).A(b^1,c)"));
    BOOST_CHECK_EQUAL(retval3[0].products()[1], Species("B(l,r^1).B(l,r^1)"));

    ReactionRule rr4;
    rr4.add_reactant(Species("A(b^1,c).A(b^1,c)"));
    rr4.add_reactant(Species("B(l,r^1).B(l,r^1)"));
    rr4.add_product(Species("A(b^1,c^2).A(b^1,c^3).B(l^2,r^4).B(l^3,r^4)"));
    rr4.set_k(1.0);

    ReactionRule::reactant_container_type reactants3(2);
    reactants3[0] = Species("A(b^1,c).A(b^1,c)");
    reactants3[1] = Species("B(l,r^1).B(l,r^1)");
    std::vector<ReactionRule> retval4 = rr4.generate(reactants3);

    BOOST_CHECK_EQUAL(retval4.size(), 2);
    BOOST_CHECK_EQUAL(retval4[0].k(), 1.0);
    BOOST_CHECK_EQUAL(retval4[1].k(), 1.0);
    BOOST_CHECK_EQUAL(retval4[0].reactants().size(), 2);
    BOOST_CHECK_EQUAL(retval4[0].reactants()[0], reactants3[0]);
    BOOST_CHECK_EQUAL(retval4[0].reactants()[1], reactants3[1]);
    BOOST_CHECK_EQUAL(retval4[1].reactants().size(), 2);
    BOOST_CHECK_EQUAL(retval4[1].reactants()[0], reactants3[0]);
    BOOST_CHECK_EQUAL(retval4[1].reactants()[1], reactants3[1]);
    BOOST_CHECK_EQUAL(retval4[0].products().size(), 1);
    BOOST_CHECK_EQUAL(retval4[1].products().size(), 1);
    BOOST_CHECK_EQUAL(format_species(retval4[0].products()[0]), Species("A(b^1,c^2).A(b^1,c^3).B(l^3,r^4).B(l^2,r^4)"));
    BOOST_CHECK_EQUAL(format_species(retval4[1].products()[0]), Species("A(b^1,c^2).A(b^1,c^3).B(l^3,r^4).B(l^2,r^4)"));

    ReactionRule rr5;
    rr5.add_reactant(Species("_(b^1)._(b^1)"));
    rr5.add_product(Species("_(b)"));
    rr5.add_product(Species("_(b)"));
    rr5.set_k(1.0);

    ReactionRule::reactant_container_type reactants4(1, Species("A(b^1).A(b^1)"));
    std::vector<ReactionRule> retval5 = rr5.generate(reactants4);

    BOOST_CHECK_EQUAL(retval5.size(), 1);
    BOOST_CHECK_EQUAL(retval5[0].k(), 1.0);
    BOOST_CHECK_EQUAL(retval5[0].reactants().size(), 1);
    BOOST_CHECK_EQUAL(retval5[0].reactants()[0], reactants4[0]);
    BOOST_CHECK_EQUAL(retval5[0].products().size(), 2);
    BOOST_CHECK_EQUAL(retval5[0].products()[0], Species("A(b)"));
    BOOST_CHECK_EQUAL(retval5[0].products()[1], Species("A(b)"));
}

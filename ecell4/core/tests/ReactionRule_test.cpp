#define BOOST_TEST_MODULE "ReactionRule_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/floating_point_comparison.hpp>

#include <ecell4/core/Species.hpp>
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

#define BOOST_TEST_MODULE "Species_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>

#include <ecell4/core/Species.hpp>
#include <ecell4/core/Context.hpp>

using namespace ecell4;


BOOST_AUTO_TEST_CASE(Species_test_constructor)
{
    Species sp1("A");
    Species sp2("A(  )");
//     Species sp2("A( a   , b ^ 1,   c)");
    Species sp3("A(a,b^1,c=p).B(a,b=u^1)");
}

BOOST_AUTO_TEST_CASE(Species_test_name)
{
    Species sp1("test");
    BOOST_CHECK_EQUAL(sp1.serial(), "test");
    Species sp2("test()");
    BOOST_CHECK_EQUAL(sp2.serial(), "test");
}

BOOST_AUTO_TEST_CASE(Species_test_attributes)
{
    Species species("test");
    species.set_attribute("attr1", "value1");
    species.set_attribute("attr2", "value2");
    BOOST_CHECK_EQUAL(species.get_attribute("attr1"), "value1");
    BOOST_CHECK_EQUAL(species.get_attribute("attr2"), "value2");
    species.remove_attribute("attr1");
}

BOOST_AUTO_TEST_CASE(Species_test_match1)
{
    Species sp1, sp2;

    sp1.add_unit(UnitSpecies("C"));
    sp1.add_unit(UnitSpecies("A"));
    sp1.add_unit(UnitSpecies("B"));
    BOOST_CHECK_EQUAL(sp1.serial(), "C.A.B");

    sp2.add_unit(UnitSpecies("A"));
    sp2.add_unit(UnitSpecies("C"));
    BOOST_CHECK_EQUAL(sp2.serial(), "A.C");

    // BOOST_CHECK(sp2.match(sp1));
    // BOOST_CHECK(!sp1.match(sp2));
}

BOOST_AUTO_TEST_CASE(Species_test_match2)
{
    Species sp1, sp2;

    sp1.add_unit(UnitSpecies("B"));
    sp1.add_unit(UnitSpecies("A"));
    sp1.add_unit(UnitSpecies("A"));
    BOOST_CHECK_EQUAL(sp1.serial(), "B.A.A");

    sp2.add_unit(UnitSpecies("B"));
    sp2.add_unit(UnitSpecies("A"));
    BOOST_CHECK_EQUAL(sp2.serial(), "B.A");

    // BOOST_CHECK(sp2.match(sp1));
    // BOOST_CHECK(!sp1.match(sp2));
}

BOOST_AUTO_TEST_CASE(Species_test_serialization)
{
    Species sp1("X(a^1).Y(a^3,b).X(a^2).Y(a^1,b^2).X(a^3)");

    BOOST_CHECK_EQUAL(
        serialize_species(Species("X(a^1).Y(a^3,b).X(a^2).Y(a^1,b^2).X(a^3)")),
        "X(a^1).Y(a^1,b).X(a^2).Y(a^2,b^3).X(a^3)");
    BOOST_CHECK_EQUAL(
        serialize_species(Species("X(a^1).Y(a^3,b^4).X(a^3).Z(a^4,b^5).Y(a^1,b^2).Z(a^2,b^5)")),
        "X(a^1).Y(a^1,b^2).Z(a^2,b^3).Z(a^4,b^3).Y(a^5,b^4).X(a^5)");
    BOOST_CHECK_EQUAL(
        serialize_species(Species("X(a^3,b^1).X(a^2,b).X(a,b^3).X(a^1,b^4).X(a^4,b^2)")),
        "X(a,b^1).X(a^1,b^2).X(a^2,b^3).X(a^3,b^4).X(a^4,b)");
}

BOOST_AUTO_TEST_CASE(Species_test_match3)
{
    BOOST_CHECK_EQUAL(count_spmatches(Species("A"), Species("A.A.A")), 3);

    BOOST_CHECK_EQUAL(count_spmatches(Species("_1._2"), Species("A.B.C")), 6);

    MatchObject::context_type::variable_container_type globals;
    globals["_1"] = "A";
    BOOST_CHECK_EQUAL(
        count_spmatches(Species("_1._2"), Species("A.B.C"), globals), 2);

    ReactionRule rr1;
    rr1.add_reactant(Species("_1(b)"));
    rr1.add_reactant(Species("_1(b)"));
    rr1.add_product(Species("_1(b^1)._1(b^1)"));

    ReactionRule::reactant_container_type reactants1;
    reactants1.push_back(Species("A(a^1,b).B(a^1,b)"));
    reactants1.push_back(Species("B(a,b)"));

    BOOST_CHECK(rrmatch(rr1, reactants1));
    BOOST_CHECK_EQUAL(count_rrmatches(rr1, reactants1), 1);

    std::vector<std::vector<Species> > retval;
    retval = rrgenerate(rr1, reactants1);
    BOOST_CHECK_EQUAL(retval.size(), 1);
    BOOST_CHECK_EQUAL(retval[0].size(), 1);
    BOOST_CHECK_EQUAL(retval[0][0].num_units(), 3);

    ReactionRule rr2;
    rr2.add_reactant(Species("A(b)"));
    rr2.add_reactant(Species("B(b)"));
    rr2.add_product(Species("A(b^1).A(b^1)"));
    ReactionRule::reactant_container_type reactants2;
    reactants2.push_back(Species("A(a^1,b).A(a^1,b)"));
    reactants2.push_back(Species("B(a^1,b).B(a^1,b^2).B(a^2,b)"));

    BOOST_CHECK(rrmatch(rr2, reactants2));
    BOOST_CHECK_EQUAL(count_rrmatches(rr2, reactants2), 4);

    ReactionRule rr3;
    rr3.add_reactant(Species("A"));
    rr3.add_product(Species("B"));
    ReactionRule::reactant_container_type reactants3;
    reactants3.push_back(Species("A"));
    retval = rrgenerate(rr3, reactants3);
    BOOST_CHECK_EQUAL(retval.size(), 1);
    BOOST_CHECK_EQUAL(retval[0].size(), 1);
    BOOST_CHECK_EQUAL(retval[0][0].serial(), "B");

    ReactionRule rr4;
    rr4.add_reactant(Species("A(b^1).B(b^1)"));
    rr4.add_product(Species("A(b)"));
    rr4.add_product(Species("B(b)"));
    ReactionRule::reactant_container_type reactants4;
    reactants4.push_back(
        Species("A(a^1,b^5).A(a^1,b^4).B(a^2,b).B(a^2,b^3).B(a^3,b^4).B(a,b^5)"));
    retval = rrgenerate(rr4, reactants4);
    BOOST_CHECK_EQUAL(retval.size(), 2);
    BOOST_CHECK_EQUAL(retval[0].size(), 2);
    BOOST_CHECK_EQUAL(retval[1].size(), 2);
    BOOST_CHECK_EQUAL(retval[0][0].num_units() + retval[0][1].num_units(), 6);
    BOOST_CHECK_EQUAL(retval[1][0].num_units() + retval[1][1].num_units(), 6);
}

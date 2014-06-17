#define BOOST_TEST_MODULE "Species_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>

#include "../Species.hpp"
#include "../Context.hpp"

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

    BOOST_CHECK(sp2.match(sp1));
    BOOST_CHECK(!sp1.match(sp2));
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

    BOOST_CHECK(sp2.match(sp1));
    BOOST_CHECK(!sp1.match(sp2));
}

// BOOST_AUTO_TEST_CASE(Species_test_get_unit)
// {
//     Species sp1;
// 
//     sp1.add_unit(UnitSpecies("B"));
//     sp1.add_unit(UnitSpecies("A"));
//     sp1.add_unit(UnitSpecies("A"));
//     BOOST_CHECK_EQUAL(sp1.serial(), "A.A.B");
// //     BOOST_CHECK_EQUAL(sp1.serial(), "B.A.A");
// 
//     BOOST_CHECK_EQUAL(sp1.get_unit(UnitSpecies("B")), 2);
//     BOOST_CHECK_EQUAL(sp1.list_sites().size(), 2);
// //     BOOST_CHECK_EQUAL(sp1.get_unit(UnitSpecies("C")), 3);
// }

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

    boost::array<Species, 2> a1 = {{Species("_1"), Species("_1")}};
    boost::array<Species, 2> b1 = {{Species("A"), Species("A.B")}};
    boost::array<Species, 2> b2 = {{Species("A"), Species("B.C")}};
    BOOST_CHECK(rrmatch(a1, b1));
    BOOST_CHECK(!rrmatch(a1, b2));
}

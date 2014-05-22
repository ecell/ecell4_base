#define BOOST_TEST_MODULE "Species_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>

#include "../Species.hpp"

using namespace ecell4;


BOOST_AUTO_TEST_CASE(Species_test_constructor)
{
    Species species("test");
}

BOOST_AUTO_TEST_CASE(Species_test_name)
{
    Species species("test");
    BOOST_CHECK_EQUAL(species.serial(), "test");
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
    // BOOST_CHECK_EQUAL(sp1.name(), "C.A.B");
    BOOST_CHECK_EQUAL(sp1.serial(), "A.B.C");
    BOOST_CHECK_EQUAL(sp1.serial(), "A.B.C");

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
    BOOST_CHECK_EQUAL(sp1.serial(), "A.A.B");

    sp2.add_unit(UnitSpecies("B"));
    sp2.add_unit(UnitSpecies("A"));
    BOOST_CHECK_EQUAL(sp2.serial(), "A.B");

    BOOST_CHECK(sp2.match(sp1));
    BOOST_CHECK(!sp1.match(sp2));
}

BOOST_AUTO_TEST_CASE(Species_test_get_unit)
{
    Species sp1;

    sp1.add_unit(UnitSpecies("B"));
    sp1.add_unit(UnitSpecies("A"));
    sp1.add_unit(UnitSpecies("A"));
    BOOST_CHECK_EQUAL(sp1.serial(), "A.A.B");
//     BOOST_CHECK_EQUAL(sp1.serial(), "B.A.A");

    BOOST_CHECK_EQUAL(sp1.get_unit(UnitSpecies("B")), 2);
    BOOST_CHECK_EQUAL(sp1.list_sites().size(), 2);
//     BOOST_CHECK_EQUAL(sp1.get_unit(UnitSpecies("C")), 3);
}

BOOST_AUTO_TEST_CASE(Species_test_serialization)
{
    Species sp1("X(a^1).Y(a^3,b).X(a^2).Y(a^1,b^2).X(a^3)");

    const std::string retval1("X(a^1).X(a^2).X(a^3).Y(a^1,b).Y(a^2,b^3)");
    // BOOST_CHECK_EQUAL(sp1.serial(), retval1);
    BOOST_CHECK_EQUAL(serialize_species(sp1), retval1);

    BOOST_CHECK_EQUAL(
        serialize_species(Species("X(a^1).Y(a^3,b^4).X(a^3).Z(a^4,b^5).Y(a^1,b^2).Z(a^2,b^5)")),
        "X(a^1).X(a^2).Y(a^1,b^3).Y(a^2,b^4).Z(a^3,b^5).Z(a^4,b^5)");
}

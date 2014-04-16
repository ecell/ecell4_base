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
    BOOST_CHECK_EQUAL(species.name(), "test");
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

BOOST_AUTO_TEST_CASE(Species_test_match)
{
    Species sp1, sp2;

    sp1.add_unit(UnitSpecies("C"));
    sp1.add_unit(UnitSpecies("A"));
    sp1.add_unit(UnitSpecies("B"));
    BOOST_CHECK_EQUAL(sp1.name(), "C.A.B");
    BOOST_CHECK_EQUAL(sp1.serial(), "A.B.C");

    sp2.add_unit(UnitSpecies("A"));
    sp2.add_unit(UnitSpecies("C"));
    BOOST_CHECK_EQUAL(sp2.name(), "A.C");
    BOOST_CHECK_EQUAL(sp2.serial(), "A.C");

    BOOST_CHECK(sp2.match(sp1));
    BOOST_CHECK(!sp1.match(sp2));
}

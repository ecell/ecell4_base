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
}

#define BOOST_TEST_MODULE "NetworkModel_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>

#include "../types.hpp"
#include "../NetworkModel.hpp"

using namespace ecell4;

BOOST_AUTO_TEST_CASE(NetworkModel_test_constructor)
{
    NetworkModel model();
}

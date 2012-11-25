#define BOOST_TEST_MODULE "hoge"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>

#include "Position3.hpp"
#include "linear_algebra.hpp"

using namespace ecell4;

BOOST_AUTO_TEST_CASE( Position3Test )
{
  Position3 pos1(1,2,3);
  BOOST_CHECK_EQUAL( pos1 * 2, Position3(2,4,6) );
}


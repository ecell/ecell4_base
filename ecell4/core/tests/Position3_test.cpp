#define BOOST_TEST_MODULE "Position3_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>

#include <ecell4/core/Position3.hpp>
#include <ecell4/core/linear_algebra.hpp>

using namespace ecell4;


BOOST_AUTO_TEST_CASE(Position3_test_multiply)
{
  Position3 pos1(1,2,3);
  BOOST_CHECK_EQUAL(pos1 * 2, Position3(2,4,6));
}

BOOST_AUTO_TEST_CASE(Position3_test_add)
{
  Position3 pos2(1,2,3);
  Position3 pos3(2,4,6);
  BOOST_CHECK_EQUAL(pos2 + pos3, Position3(3,6,9));
}

BOOST_AUTO_TEST_CASE(Position3_test_sub)
{
  Position3 pos4(2,4,6);
  Position3 pos5(1,2,3);
  BOOST_CHECK_EQUAL(pos4 - pos5, Position3(1,2,3));
}

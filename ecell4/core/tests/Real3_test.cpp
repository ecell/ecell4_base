#define BOOST_TEST_MODULE "Real3_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <ecell4/core/Real3.hpp>
#include <ecell4/core/linear_algebra.hpp>

using namespace ecell4;


BOOST_AUTO_TEST_CASE(Real3_test_multiply)
{
  Real3 pos1(1,2,3);
  BOOST_CHECK_EQUAL(pos1 * 2, Real3(2,4,6));
}

BOOST_AUTO_TEST_CASE(Real3_test_add)
{
  Real3 pos2(1,2,3);
  Real3 pos3(2,4,6);
  BOOST_CHECK_EQUAL(pos2 + pos3, Real3(3,6,9));
}

BOOST_AUTO_TEST_CASE(Real3_test_sub)
{
  Real3 pos4(2,4,6);
  Real3 pos5(1,2,3);
  BOOST_CHECK_EQUAL(pos4 - pos5, Real3(1,2,3));
}

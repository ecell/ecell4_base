#define BOOST_TEST_MODULE "array_helper_test"

#include <boost/test/included/unit_test.hpp>
#include "array_helper.hpp"

BOOST_AUTO_TEST_CASE(basic)
{
    boost::array<int, 3> a(array_gen<int>(1, 2, 3));
    BOOST_CHECK_EQUAL(1, a[0]);
    BOOST_CHECK_EQUAL(2, a[1]);
    BOOST_CHECK_EQUAL(3, a[2]);
}

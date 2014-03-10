#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#define BOOST_TEST_MODULE "Vector3_test"

#include <boost/test/included/unit_test.hpp>
#include "Vector3.hpp"
#include "geometry.hpp"

BOOST_AUTO_TEST_CASE(test_cyclic_transpose)
{
    typedef Vector3<double> vec;

    BOOST_CHECK_EQUAL(
        cyclic_transpose(vec(8., 8., 8.), vec(1., 1., 1.), 10.),
        vec(-2., -2., -2.));
}

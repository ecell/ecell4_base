#define BOOST_TEST_MODULE "linear_algebra_test"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "utils/array_helper.hpp"
#include "geometry.hpp"

BOOST_AUTO_TEST_CASE(test_rotate_vector)
{
    {
        boost::array<double, 3> const r(rotate_vector(array_gen(1., 0., 0.),
                                              array_gen(1., 0., 0.), .1));
        BOOST_CHECK_EQUAL(r[0], 1.);
        BOOST_CHECK_EQUAL(r[1], 0.);
        BOOST_CHECK_EQUAL(r[2], 0.);
    }

    {
        boost::array<double, 3> const r(rotate_vector(array_gen(1., 0., 0.),
                                              array_gen(0., 0., 1.), M_PI));
        BOOST_CHECK_CLOSE(r[0], -1., 0.000000001);
        // r[1]{1.2246063538223773e-16}
        BOOST_CHECK_SMALL(r[1], 0.000000001);
        BOOST_CHECK_CLOSE(r[2], 0., 0.0000000001);
    }
}

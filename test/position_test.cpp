#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#define BOOST_TEST_MODULE "position"

#include <boost/config.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/mpl/list.hpp>

#include "position.hpp"

typedef boost::mpl::list3<long double, double, float> scalar_types;

template<typename T_>
struct get_tolerance
{
};

template<>
struct get_tolerance<float>
{
    BOOST_STATIC_CONSTANT(float, value = 1e-8);
};

template<>
struct get_tolerance<double>
{
    BOOST_STATIC_CONSTANT(double, value = 1e-17);
};

template<>
struct get_tolerance<long double>
{
    BOOST_STATIC_CONSTANT(long double, value = 1e-26);
};

BOOST_AUTO_TEST_CASE_TEMPLATE(distance_sq, T_, scalar_types)
{
    typedef ::position<T_> position;

    BOOST_CHECK_CLOSE((T_)3.0f,
            position(1, 1, 1).distance_sq(position(0, 0, 0)),
            get_tolerance<T_>::value);
    BOOST_CHECK_CLOSE((T_)3.0f,
            position(1, 2, 1).distance_sq(position(0, 1, 0)),
            get_tolerance<T_>::value);
    BOOST_CHECK_CLOSE((T_)0.75f,
            position(0.5, 0.5, 1).distance_sq(position(0, 1, 0.5)),
            get_tolerance<T_>::value);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(add_op, T_, scalar_types)
{
    typedef ::position<T_> position;

    BOOST_CHECK_EQUAL(position(1, 1, 1),
            position(1, 1, 1) + position(0, 0, 0));
    BOOST_CHECK_EQUAL(position(1, 3, 1),
            position(1, 2, 1) + position(0, 1, 0));
    BOOST_CHECK_EQUAL(position(0.5, 1.5, 1.5),
            position(0.5, 0.5, 1) + position(0, 1, 0.5));
}


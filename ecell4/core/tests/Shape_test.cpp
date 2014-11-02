#define BOOST_TEST_MODULE "Shape_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/floating_point_comparison.hpp>

#include <ecell4/core/Sphere.hpp>

using namespace ecell4;

struct Fixture
{
    const Real3 center;
    const Real radius;
    Sphere sphere;
    Fixture() :
        center(2.5e-6, 2.5e-6, 2.5e-6),
        radius(2.5e-7), sphere(center, radius)
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(suite, Fixture)

BOOST_AUTO_TEST_CASE(Shape_test_constructor)
{
}

BOOST_AUTO_TEST_CASE(Shape_test_is_inside)
{
    BOOST_CHECK(sphere.is_inside(center) <= 0);
    BOOST_CHECK(sphere.is_inside(Real3(2.3e-6, 2.5e-6, 2.5e-6)) <= 0);
    BOOST_CHECK(sphere.is_inside(Real3(2.5e-6, 2.3e-6, 2.5e-6)) <= 0);
    BOOST_CHECK(sphere.is_inside(Real3(2.5e-6, 2.5e-6, 2.3e-6)) <= 0);
    BOOST_CHECK(sphere.is_inside(Real3(2.2e-6, 2.5e-6, 2.5e-6)) > 0);
    BOOST_CHECK(sphere.is_inside(Real3(2.5e-6, 2.2e-6, 2.5e-6)) > 0);
    BOOST_CHECK(sphere.is_inside(Real3(2.5e-6, 2.5e-6, 2.2e-6)) > 0);
}

BOOST_AUTO_TEST_SUITE_END()

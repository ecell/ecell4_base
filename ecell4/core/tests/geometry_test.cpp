#define BOOST_TEST_MODULE "geometry_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/geometry.hpp>

using namespace ecell4;

static const Real tolerance = 1.e-10;
static const std::size_t  N = 10000;

BOOST_AUTO_TEST_CASE(geometry_test_rotate)
{
    std::shared_ptr<RandomNumberGenerator> rng(new GSLRandomNumberGenerator());
    const Real3 unitx(1, 0, 0);
    const Real3 unity(0, 1, 0);
    const Real3 unitz(0, 0, 1);

    for(std::size_t i=0; i<N; ++i)
    {
        const Real angle = rng->uniform(-M_PI, M_PI);
        const Real3 result = rotate(angle, unitz, unitx);

        BOOST_CHECK_CLOSE_FRACTION(result[0], std::cos(angle), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(result[1], std::sin(angle), tolerance);
        BOOST_CHECK_SMALL(result[2], tolerance);
    }
}

BOOST_AUTO_TEST_CASE(geometry_test_angle)
{
    std::shared_ptr<RandomNumberGenerator> rng(new GSLRandomNumberGenerator());
    const Real3 unitz(0, 0, 1);
    const Real3 unitx(1, 0, 0);

    for(std::size_t i=0; i<N; ++i)
    {
        const Real angle = rng->uniform(-M_PI, M_PI);
        const Real a = rng->uniform(0, 2);
        const Real b = rng->uniform(0, 2);
        const Real3 rotated = rotate(angle, unitz, unitx);
        const Real result = angle(unitx * a, rotated * b);

        BOOST_CHECK_CLOSE_FRACTION(result, angle, tolerance);
    }
}


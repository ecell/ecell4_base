#define BOOST_TEST_MODULE "Barycentric_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Barycentric.hpp>

using namespace ecell4;

static const Real tolerance = 1.e-10;
static const std::size_t  N = 10000;

BOOST_AUTO_TEST_CASE(Barycentric_test_operator)
{
    const Barycentric b1(1., 2., 3.);
    const Barycentric b2(6., 4., 5.);
    const Barycentric b3 = b1 + b2;

    BOOST_CHECK_CLOSE_FRACTION(b3[0], 7.0, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(b3[1], 6.0, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(b3[2], 8.0, tolerance);

    const Barycentric b4 = b2 - b1;

    BOOST_CHECK_CLOSE_FRACTION(b4[0], 5.0, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(b4[1], 2.0, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(b4[2], 2.0, tolerance);
}


BOOST_AUTO_TEST_CASE(Barycentric_test_is_inside)
{
    const Barycentric b1(1., 0., 0.);
    const Barycentric b2(0.2, 0.3, 0.5);

    BOOST_CHECK(is_inside(b1));
    BOOST_CHECK(is_inside(b2));

    const Barycentric b3(0., 0., 0.);
    const Barycentric b4(1., 2., 3.);
    const Barycentric b5(1., 1., -1.);

    BOOST_CHECK(!is_inside(b3));
    BOOST_CHECK(!is_inside(b4));
    BOOST_CHECK(!is_inside(b5));
}

BOOST_AUTO_TEST_CASE(Barycentric_test_on_plane)
{
    const Barycentric b1(1., 0., 0.);
    const Barycentric b2(0.2, 0.3, 0.5);
    const Barycentric b3(1., 1., -1.);

    BOOST_CHECK(on_plane(b1));
    BOOST_CHECK(on_plane(b2));
    BOOST_CHECK(on_plane(b3));

    const Barycentric b4(0., 0., 0.);
    const Barycentric b5(1., 2., 3.);

    BOOST_CHECK(!on_plane(b4));
    BOOST_CHECK(!on_plane(b5));
}

BOOST_AUTO_TEST_CASE(Barycentric_test_transformation)
{
    std::shared_ptr<RandomNumberGenerator> rng(new GSLRandomNumberGenerator());

    for(std::size_t i=0; i<N; ++i)
    {
        const Real a = rng->uniform(0., 1.);
        const Real b = rng->uniform(0., 1. - a);
        const Real c = 1. - a - b;
        const Barycentric bary(a, b, c);

        const Real3 v0(0., 0., 0.);
        const Real3 v1 = rng->direction3d(1.);
        Real3 v2;
        while(true)
        {
            v2 = rng->direction3d(1.);
            const Real dot = std::abs(dot_product(v1, v2));
            if(std::abs(dot - 1.0) > tolerance) break;
        }
        const Triangle tri(v0, v1, v2);

        const Real3 absolute = to_absolute(bary, tri);
        const Barycentric ret = to_barycentric(absolute, tri);

        BOOST_CHECK_CLOSE_FRACTION(bary[0], ret[0], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(bary[1], ret[1], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(bary[2], ret[2], tolerance);
    }
}


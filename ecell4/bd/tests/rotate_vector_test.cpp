#define BOOST_TEST_MODULE "rotate_vector_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/bd/rotate_vector.hpp>

using namespace ecell4;

static const Real tolerance = 1.e-10;
static const std::size_t  N = 10000;

BOOST_AUTO_TEST_CASE(rotate_test)
{
    const Real3 x(1., 0., 0.);
    const Real3 n(0., 0., 1.);

    const Real dtheta = M_PI * 2. / N;
    for(std::size_t i=0; i<N; ++i)
    {
        const Real theta = i * dtheta;
        const Real3 vec(rotate(theta, n, x));

        if(std::abs(std::cos(theta)) < tolerance)
            BOOST_CHECK_SMALL(vec[0], tolerance);
        else
            BOOST_CHECK_CLOSE_FRACTION(vec[0], std::cos(theta), tolerance);

        if(std::abs(std::sin(theta)) < tolerance)
            BOOST_CHECK_SMALL(vec[1], tolerance);
        else
            BOOST_CHECK_CLOSE_FRACTION(vec[1], std::sin(theta), tolerance);

        BOOST_CHECK_SMALL(vec[2], tolerance);
    }
}

BOOST_AUTO_TEST_CASE(angle_test)
{
    const Real3 x(1., 0., 0.);
    const Real3 n(0., 0., 1.);

    const std::size_t Nh = N / 2;
    const Real dtheta = M_PI / Nh;
    for(std::size_t i=0; i<Nh; ++i)
    {
        const Real theta = i * dtheta;
        const Real3 vec(rotate(theta, n, x));

        const Real phi = calc_angle(x, vec);
        BOOST_CHECK_CLOSE_FRACTION(phi, theta, 1e-8);
    }

    for(std::size_t i=Nh+1; i<N; ++i)
    {
        const Real theta = i * dtheta;
        const Real3 vec(rotate(theta, n, x));

        const Real phi = calc_angle(x, vec);
        BOOST_CHECK_CLOSE_FRACTION(phi, 2 * M_PI - theta, 1e-8);
    }
}



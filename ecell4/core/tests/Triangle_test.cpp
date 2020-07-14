#define BOOST_TEST_MODULE "Triangle_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/BoundaryCondition.hpp>
#include <ecell4/core/Triangle.hpp>
#include <random>

using ecell4::Real;
using ecell4::Real3;
using ecell4::Triangle;

BOOST_AUTO_TEST_CASE(Triangle_distance_test)
{
    constexpr std::size_t N = 100;
    constexpr Real tol = 1e-8;
    // test triangle
    //
    //     y   .'
    //  A  | .'
    // ----|'     B
    //    1|`.         .'
    //  F  |  `.     .'
    //     | G  `. .'
    // ----+------`----- x
    //    0|     1|
    //  E  |   D  |  C

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<Real> uni(0.0, 1.0);

    ecell4::UnlimitedBoundary boundary;

    const Triangle tri(Real3(0, 0, 0), Real3(1, 0, 0), Real3(0, 1, 0));
    BOOST_TEST_MESSAGE("triangle constructed");

    // check region A
    for(std::size_t i=0; i<N; ++i)
    {
        const auto x = -uni(mt);
        const auto y = (1.0 - x) + uni(mt);
        const auto z = 2 * (uni(mt) - 0.5);
        const Real3 pos(x,y,z);

        const auto dist = distance_point_Triangle(pos, tri, boundary);
        BOOST_CHECK_CLOSE_FRACTION(dist, length(pos - Real3(0, 1, 0)), tol);
    }

    // check region B
    for(std::size_t i=0; i<N; ++i)
    {
        const auto x = uni(mt);
        const auto y = (1.0 - x) + uni(mt) * 2 * x;
        const auto z = 2 * (uni(mt) - 0.5);
        const Real3 pos(x,y,z);

        const auto ref_dist = std::sqrt(0.5 * (x+y-1) * (x+y-1) + z * z);

        const auto dist = distance_point_Triangle(pos, tri, boundary);
        BOOST_CHECK_CLOSE_FRACTION(dist, ref_dist, tol);
    }

    // check region C
    for(std::size_t i=0; i<N; ++i)
    {
        const auto y = -uni(mt);
        const auto x = (1.0 - y) + uni(mt);
        const auto z = 2 * (uni(mt) - 0.5);
        const Real3 pos(x,y,z);

        const auto dist = distance_point_Triangle(pos, tri, boundary);
        BOOST_CHECK_CLOSE_FRACTION(dist, length(pos - Real3(1, 0, 0)), tol);
    }

    // check region D
    for(std::size_t i=0; i<N; ++i)
    {
        const auto x = uni(mt);
        const auto y = -uni(mt);
        const auto z = 2 * (uni(mt) - 0.5);
        const Real3 pos(x,y,z);

        const auto dist = distance_point_Triangle(pos, tri, boundary);
        BOOST_CHECK_CLOSE_FRACTION(dist, std::sqrt(y*y + z*z), tol);
    }

    // check region E
    for(std::size_t i=0; i<N; ++i)
    {
        const auto x = -uni(mt);
        const auto y = -uni(mt);
        const auto z = 2 * (uni(mt) - 0.5);
        const Real3 pos(x,y,z);

        const auto dist = distance_point_Triangle(pos, tri, boundary);
        BOOST_CHECK_CLOSE_FRACTION(dist, length(pos - Real3(0, 0, 0)), tol);
    }

    // check region F
    for(std::size_t i=0; i<N; ++i)
    {
        const auto x = -uni(mt);
        const auto y = uni(mt);
        const auto z = 2 * (uni(mt) - 0.5);
        const Real3 pos(x,y,z);

        const auto dist = distance_point_Triangle(pos, tri, boundary);
        BOOST_CHECK_CLOSE_FRACTION(dist, std::sqrt(x*x + z*z), tol);
    }

    // check region G
    for(std::size_t i=0; i<N; ++i)
    {
        const auto x = uni(mt);
        const auto y = uni(mt) * (1.0 - x);
        const auto z = 2 * (uni(mt) - 0.5);
        const Real3 pos(x,y,z);

        const auto dist = distance_point_Triangle(pos, tri, boundary);
        BOOST_CHECK_CLOSE_FRACTION(dist, std::abs(z), tol);
    }
}


BOOST_AUTO_TEST_CASE(Triangle_distance_under_PBC)
{
    constexpr std::size_t N = 100;
    constexpr Real tol = 1e-8;

    // A nasty test case where the point which is the nearest periodic image to
    // the center of mass of the triangle is not the nearest periodic image to
    // the triangle itself.

    const Real3 edge_lengths(1.0, 20.0, 20.0);

    const ecell4::PeriodicBoundary  periodic(edge_lengths);

    const Real3 probe(0.1, 10.0, 10.0);

    const Triangle tri(
            Real3(0.1 + 3.0 / 4.0, 10.0 +  0.0,     10.0 + -std::sqrt(3.0) / 12.0),
            Real3(0.1 + 1.0 / 4.0, 10.0 + -1.0/3.0, 10.0 + -std::sqrt(3.0) /  4.0),
            Real3(0.1 + 1.0 / 4.0, 10.0 +  1.0/3.0, 10.0 + -std::sqrt(3.0) /  4.0));

    const Real3 CoM = (tri.vertex_at(0) + tri.vertex_at(1) + tri.vertex_at(2)) / 3.0;

    const Real3 nearest_to_CoM = periodic.periodic_transpose(probe, CoM);

    // check that probe is the nearest image to the CoM
    BOOST_CHECK_CLOSE_FRACTION(nearest_to_CoM[0], probe[0], tol);
    BOOST_CHECK_CLOSE_FRACTION(nearest_to_CoM[1], probe[1], tol);
    BOOST_CHECK_CLOSE_FRACTION(nearest_to_CoM[2], probe[2], tol);

    const auto dist_noPBC = distance_point_Triangle(probe, tri);
    const auto dist_PBC   = distance_point_Triangle(probe, tri, periodic);

    // check if distance_point_Triangle finds the nearest periodic image
    BOOST_CHECK_CLOSE_FRACTION(dist_noPBC, 0.5, tol);
    BOOST_CHECK_CLOSE_FRACTION(dist_PBC, std::sqrt(3.0) / 6.0, tol);

    BOOST_CHECK(dist_PBC < dist_noPBC);
}

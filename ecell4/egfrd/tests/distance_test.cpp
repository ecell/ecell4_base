#define BOOST_TEST_MODULE "distance_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include "../TriangleOperation.hpp"
#include "../BarycentricCoordinate.hpp"
#include "../FaceTriangle.hpp"
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Real3.hpp>
#include <cmath>

BOOST_AUTO_TEST_CASE(closest_point_is_inside_of_triangle_test)
{
    std::size_t test_times = 1000;

    ecell4::GSLRandomNumberGenerator rng;

    for(std::size_t i=0; i<test_times; ++i)
    {
        const ecell4::Real3 center(rng.random(), rng.random(), rng.random());
        const ecell4::Real theta = rng.uniform(0., M_PI);
        const ecell4::Real phi = rng.uniform(0., 2 * M_PI);
        const ecell4::Real3 direction(
            sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
        BOOST_CHECK_CLOSE_FRACTION(length(direction), 1e0, 1e-12);

        const ecell4::Real distance = std::abs(rng.random());
        const ecell4::Real3 closest_position = center + direction * distance;

        ecell4::Real3 random(rng.random(), rng.random(), rng.random());
        ecell4::Real3 on_plane_ = ecell4::cross_product(direction, center);
        while(length(on_plane_) == 0.)
        {
            random = ecell4::Real3(rng.random(), rng.random(), rng.random());
            on_plane_ = ecell4::cross_product(direction, center);
        }

        const ecell4::Real3 on_plane = on_plane_;
        const ecell4::Real dot_a = ecell4::dot_product(on_plane, direction);
        BOOST_CHECK_SMALL(dot_a, 1e-12);

        ecell4::Real alpha_ = rng.uniform(0., M_PI);
        while(alpha_ == 0.){alpha_ = rng.uniform(0., M_PI);}
        const ecell4::Real alpha = alpha_;

        ecell4::Real beta_ = rng.uniform(M_PI, alpha + M_PI);
        while(alpha == beta_ || beta_ == 2 * M_PI){beta_ = rng.uniform(M_PI, alpha + M_PI);}
        const ecell4::Real beta = beta_;

        const ecell4::Real3 direction_b = rotation(alpha, direction, on_plane);
        const ecell4::Real dot_b = ecell4::dot_product(direction_b, direction);
        BOOST_CHECK_SMALL(dot_b, 1e-12);

        const ecell4::Real3 direction_c = rotation(beta, direction, on_plane);
        const ecell4::Real dot_c = ecell4::dot_product(direction_c, direction);
        BOOST_CHECK_SMALL(dot_c, 1e-12);

        const ecell4::Real3 a = closest_position + on_plane * rng.random();
        const ecell4::Real3 b = closest_position + direction_b * rng.random();
        const ecell4::Real3 c = closest_position + direction_c * rng.random();

        boost::array<ecell4::Real3, 3> arr;
        arr[0] = a;
        arr[1] = b;
        arr[2] = c;

        // closest_position is on the triangle
        const FaceTriangle<ecell4::Real3> tri(arr);
        const Barycentric<ecell4::Real> bary = make_barycentric(closest_position, tri);
        const ecell4::Real3 absolute = make_absolute(bary, tri);
        BOOST_CHECK_CLOSE_FRACTION(closest_position[0], absolute[0], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(closest_position[1], absolute[1], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(closest_position[2], absolute[2], 1e-12);
        // end

        const ecell4::Real3 result = ecell4::closest_point(center, arr);
        const ecell4::Real result_dist = ecell4::length(result - center);

        // result is on the triangle
        const Barycentric<ecell4::Real> result_bary = make_barycentric(result, tri);
        const ecell4::Real3 result_absolute = make_absolute(result_bary, tri);
        BOOST_CHECK_CLOSE_FRACTION(result[0], result_absolute[0], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(result[1], result_absolute[1], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(result[2], result_absolute[2], 1e-12);
        // end

        BOOST_CHECK_CLOSE_FRACTION(result_dist, distance, 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(closest_position[0], result[0], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(closest_position[1], result[1], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(closest_position[2], result[2], 1e-12);
    }
}

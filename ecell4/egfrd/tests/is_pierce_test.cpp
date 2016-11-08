#define BOOST_TEST_MODULE "is_pierce_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include "../TriangleOperation.hpp"
#include "../FaceTriangle.hpp"
#include "../BarycentricCoordinate.hpp"
#include "../geometry.hpp"
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Real3.hpp>
#include <cmath>


BOOST_AUTO_TEST_CASE(test_bare_is_pierce)
{
    std::size_t test_times = 1000;

    ecell4::GSLRandomNumberGenerator rng;

    for(std::size_t i=0; i<test_times; ++i)
    {
        const ecell4::Real3 begin(rng.random(), rng.random(), rng.random());

        const ecell4::Real theta = rng.uniform(0., M_PI);
        const ecell4::Real phi = rng.uniform(0., 2 * M_PI);
        const ecell4::Real3 direction(
            sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
        BOOST_CHECK_CLOSE(length(direction), 1e0, 1e-12);

        const ecell4::Real distance = std::abs(rng.random());
        const ecell4::Real3 through = begin + direction * distance;
        const ecell4::Real3 end = begin + direction * (2 * distance);

        ecell4::Real3 axis_(rng.random(), rng.random(), rng.random());
        while(std::abs(ecell4::dot_product(axis_, direction)) == 1 || length(axis_) == 0.)
        {
            axis_ = ecell4::Real3(rng.random(), rng.random(), rng.random());
        }
        const ecell4::Real3 axis = axis_ / length(axis_);

        ecell4::Real3 random(rng.random(), rng.random(), rng.random());
        ecell4::Real3 on_plane_ = ecell4::cross_product(axis, random);
        while(length(on_plane_) == 0.)
        {
            random = ecell4::Real3(rng.random(), rng.random(), rng.random());
            on_plane_ = ecell4::cross_product(axis, random);
        }

        const ecell4::Real3 on_plane = on_plane_;
        const ecell4::Real dot_a = ecell4::dot_product(on_plane, axis);
        BOOST_CHECK_SMALL(dot_a, 1e-12);

        ecell4::Real alpha_ = rng.uniform(0., M_PI);
        while(alpha_ == 0.){alpha_ = rng.uniform(0., M_PI);}
        const ecell4::Real alpha = alpha_;

        ecell4::Real beta_ = rng.uniform(M_PI, alpha + M_PI);
        while(alpha == beta_ || beta_ == 2 * M_PI){beta_ = rng.uniform(M_PI, alpha + M_PI);}
        const ecell4::Real beta = beta_;

        const ecell4::Real3 direction_b = rotate_vector(on_plane, axis, alpha);
        const ecell4::Real dot_b = ecell4::dot_product(direction_b, axis);
        BOOST_CHECK_SMALL(dot_b, 1e-12);

        const ecell4::Real3 direction_c = rotate_vector(on_plane, axis, beta);
        const ecell4::Real dot_c = ecell4::dot_product(direction_c, axis);
        BOOST_CHECK_SMALL(dot_c, 1e-12);

        const ecell4::Real3 a = through + on_plane * rng.random();
        const ecell4::Real3 b = through + direction_b * rng.random();
        const ecell4::Real3 c = through + direction_c * rng.random();

        boost::array<ecell4::Real3, 3> arr;
        arr[0] = a;
        arr[1] = b;
        arr[2] = c;
        boost::array<ecell4::Real3, 3> rev;
        rev[0] = a;
        rev[2] = b;
        rev[1] = c;

        // closest_position is on the triangle
        const FaceTriangle<ecell4::Real3> tri(arr);
        const Barycentric<ecell4::Real> bary = make_barycentric(through, tri);
        const ecell4::Real3 absolute = make_absolute(bary, tri);
        BOOST_CHECK_CLOSE_FRACTION(through[0], absolute[0], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(through[1], absolute[1], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(through[2], absolute[2], 1e-12);
        // end

        {
        const std::pair<bool, ecell4::Real3> result = is_pierce(begin, end, tri);
        BOOST_CHECK(result.first);
        BOOST_CHECK_CLOSE_FRACTION(result.second[0], through[0], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(result.second[1], through[1], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(result.second[2], through[2], 1e-12);
        }

        {
        const FaceTriangle<ecell4::Real3> rv(rev);
        const std::pair<bool, ecell4::Real3> result = is_pierce(begin, end, rv);
        BOOST_CHECK(result.first);
        BOOST_CHECK_CLOSE_FRACTION(result.second[0], through[0], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(result.second[1], through[1], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(result.second[2], through[2], 1e-12);
        }
    }
}

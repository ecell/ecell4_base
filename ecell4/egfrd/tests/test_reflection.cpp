#define BOOST_TEST_MODULE "reflection_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include "../Vector3Operation.hpp"
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Real3.hpp>

BOOST_AUTO_TEST_CASE(BarycentricCoordinate_test)
{
    const std::size_t test_times = 10000;
    ecell4::GSLRandomNumberGenerator rng;

    for(std::size_t i=0; i<test_times; ++i)
    {
        const ecell4::Real3 center(rng.random(), rng.random(), rng.random());
        const ecell4::Real theta = rng.uniform(0., M_PI);
        const ecell4::Real phi = rng.uniform(0., 2 * M_PI);
        const ecell4::Real3 normal(
            sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
        BOOST_CHECK_CLOSE_FRACTION(length(normal), 1e0, 1e-12);

        const ecell4::Real3 begin(rng.random(), rng.random(), rng.random());
        const ecell4::Real3 end(rng.random(), rng.random(), rng.random());

        const ecell4::Real3 reflected = reflect_plane(begin, end, normal, center);

        const ecell4::Real norm_before = ecell4::dot_product(begin - center, normal);
        const ecell4::Real norm_after  = ecell4::dot_product(reflected - center, normal);

        const bool side_before = norm_before > 0.;
        const bool side_after  = norm_after > 0.;

        if(side_before != side_after)
        {
            std::cout << "norm_bef: " << norm_before << std::endl;
            std::cout << "norm_aft: " << norm_after  << std::endl;
            std::cout << "begin   : " << begin << std::endl;
            std::cout << "end     : " << end << std::endl;
        }

        BOOST_CHECK_EQUAL(side_before, side_after);
    }
}


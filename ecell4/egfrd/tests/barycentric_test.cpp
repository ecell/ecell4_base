#define BOOST_TEST_MODULE "GreensFunction3DRadInf_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include "../BarycentricCoordinate.hpp"
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Real3.hpp>

BOOST_AUTO_TEST_CASE(BarycentricCoordinate_test)
{
    ecell4::GSLRandomNumberGenerator rng;

    for(std::size_t i=0; i<100; ++i)
    {
        ecell4::Real3 a(rng.random(), 0., 0.);
        ecell4::Real3 b(0., rng.random(), 0.);
        ecell4::Real3 c(0., 0., rng.random());
        while(a[0] == 0 && b[1] == 0){a[0] = rng.random();}
        while(b[1] == 0 && c[2] == 0){b[0] = rng.random();}
        while(c[2] == 0 && a[0] == 0){c[0] = rng.random();}
        FaceTriangle<ecell4::Real3> tri(a, b, c);

        const ecell4::Real3 pos =
            a + (b - a) * rng.random() + (c - a) * rng.random();

        const Barycentric<ecell4::Real> bary = make_barycentric(pos, tri);
        const ecell4::Real3 result = make_absolute(bary, tri);

        BOOST_CHECK_CLOSE_FRACTION(pos[0], result[0], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(pos[1], result[1], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(pos[2], result[2], 1e-12);
    }
}


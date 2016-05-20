#define BOOST_TEST_MODULE "GreensFunction3DRadInf_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include "../GreensFunction3DRadInf.hpp"

using namespace ecell4;


BOOST_AUTO_TEST_CASE(GreensFunction3DRadInf_test_constructor)
{
    const Real D = 2e-12;
    const Real kf = 0;
    const Real r0 = 1.0084e-08;
    const Real sigma = 1e-08;

    GreensFunction3DRadInf gf(D, kf, r0, sigma);
}

// BOOST_AUTO_TEST_CASE(GreensFunction3DRadInf_test_drawTheta)
// {
//     const Real D = 2e-12;
//     const Real kf = 0;
//     const Real r0 = 1.0084e-08;
//     const Real sigma = 1e-08;
// 
//     GreensFunction3DRadInf gf(D, kf, r0, sigma);
// 
//     const Real r = 1.11944e-08;
//     const Real t = 5.08006e-08;
//     const Real rnd = 0.5;  //XXX
// 
//     const Real theta = gf.drawTheta(rnd, r, t);
//     BOOST_CHECK(theta >= 0);
// }

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
// gsl_integration_qag =>
// tol = 4.45456e+13
// n = 4
// r = 1.11944e-08
// t = 5.08006e-08

//     const Real L(1e-6);
//     const Real3 input(L, L, L);
//     const Integer3 matrix_sizes(3, 3, 3);
//     boost::shared_ptr<RandomNumberGenerator> rng(new GSLRandomNumberGenerator());
// 
//     GreensFunction3DRadInf target(input, matrix_sizes, rng);
// 
//     const Real3& output(target.edge_lengths());
//     for (Real3::size_type dim(0); dim < 3; ++dim)
//     {
//         BOOST_CHECK(output[dim] > 0);
//         BOOST_CHECK_EQUAL(output[dim], input[dim]);
//     }
// }

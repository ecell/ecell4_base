#define BOOST_TEST_MODULE "boosttest"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "../Position3.hpp"
#include "../linear_algebra.hpp"
#include "../types.hpp"
#include "../CompartmentSpace.hpp"
#include "../RandomNumberGenerator.hpp"

using namespace ecell4;

// template<typename Timpl_2>
// void RandomNumberGenerator_test_seed_template()
// {
//   Timpl_2 target2;
//   target2.seed(100);
// }

// BOOST_AUTO_TEST_CASE(RandomNumberGenerator_test_seed)
// {
//   RandomNumberGenerator_test_seed_template<GSLRandomNumberGenerator>();
// }

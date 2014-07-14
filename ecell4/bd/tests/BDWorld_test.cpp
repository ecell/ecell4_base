#define BOOST_TEST_MODULE "BDWorld_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include "../BDWorld.hpp"

using namespace ecell4;
using namespace ecell4::bd;


BOOST_AUTO_TEST_CASE(BDWorld_test_constructor)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    boost::shared_ptr<RandomNumberGenerator> rng(new GSLRandomNumberGenerator());

    BDWorld target(edge_lengths, rng);
}

BOOST_AUTO_TEST_CASE(BDWorld_test_edge_lengths)
{
    const Real L(1e-6);
    const Position3 input(L, L, L);
    boost::shared_ptr<RandomNumberGenerator> rng(new GSLRandomNumberGenerator());

    BDWorld target(input, rng);

    const Position3& output(target.edge_lengths());
    for (Position3::size_type dim(0); dim < 3; ++dim)
    {
        BOOST_CHECK(output[dim] > 0);
        BOOST_CHECK_EQUAL(output[dim], input[dim]);
    }
}

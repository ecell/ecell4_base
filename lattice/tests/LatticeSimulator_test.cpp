#define BOOST_TEST_MODULE "LatticeSimulator_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "../LatticeSimulator.hpp"
#include "../../core/SerialIDGenerator.hpp"

using namespace ecell4;
using namespace ecell4::lattice;

BOOST_AUTO_TEST_CASE(LatticeSimulator_test_constructor)
{
    boost::shared_ptr<NetworkModel> model;
    boost::shared_ptr<GSLRandomNumberGenerator> rng;
    Position3 edge_lengths(1e-6, 1e-6, 1e-6);
    boost::shared_ptr<LatticeWorld> world(new LatticeWorld(edge_lengths, rng));
    LatticeSimulator sim(model, world);
}

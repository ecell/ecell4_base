#define BOOST_TEST_MODULE "BDSimulator_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>

#include <ecell4/core/NetworkModel.hpp>
#include "../BDSimulator.hpp"

using namespace ecell4;
using namespace ecell4::bd;


BOOST_AUTO_TEST_CASE(BDSimulator_test_constructor)
{
    Real const L(1e-6);
    Position3 const edge_lengths(L, L, L);

    boost::shared_ptr<Model> model(new NetworkModel());
    boost::shared_ptr<BDWorld> world(new BDWorld(edge_lengths));
    GSLRandomNumberGenerator rng;

    BDSimulator target(model, world, rng);
}

BOOST_AUTO_TEST_CASE(BDSimulator_test_step)
{
    Real const L(1e-6);
    Position3 const edge_lengths(L, L, L);

    boost::shared_ptr<Model> model(new NetworkModel());
    boost::shared_ptr<BDWorld> world(new BDWorld(edge_lengths));
    GSLRandomNumberGenerator rng;

    BDSimulator target(model, world, rng);
    target.step();
}

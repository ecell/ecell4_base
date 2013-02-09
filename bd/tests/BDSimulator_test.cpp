#define BOOST_TEST_MODULE "BDSimulator_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>

#include <ecell4/core/NetworkModel.hpp>
#include "../BDSimulator.hpp"

using namespace ecell4;
using namespace ecell4::bd;


BOOST_AUTO_TEST_CASE(BDSimulator_test_constructor)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    boost::shared_ptr<RandomNumberGenerator> rng(new GSLRandomNumberGenerator());

    boost::shared_ptr<Model> model(new NetworkModel());
    boost::shared_ptr<BDWorld> world(new BDWorld(edge_lengths, rng));

    BDSimulator target(model, world);
}

BOOST_AUTO_TEST_CASE(BDSimulator_test_step)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    boost::shared_ptr<RandomNumberGenerator> rng(new GSLRandomNumberGenerator());

    boost::shared_ptr<Model> model(new NetworkModel());
    boost::shared_ptr<BDWorld> world(new BDWorld(edge_lengths, rng));

    BDSimulator target(model, world);
    target.step();
}

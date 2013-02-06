#define BOOST_TEST_MODULE "EGFRDWorld_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>

#include <ecell4/core/NetworkModel.hpp>

#include "../EGFRDWorld.hpp"
#include "../EGFRDSimulatorWrapper.hpp"

using namespace ecell4;
using namespace ecell4::egfrd;


BOOST_AUTO_TEST_CASE(EGFRDWorld_test_constructor)
{
    boost::shared_ptr<NetworkModel> model(new NetworkModel());

    Real const world_size(1e-6);
    boost::shared_ptr<EGFRDWorld> world(new EGFRDWorld(world_size));

    EGFRDSimulatorWrapper sim(model, world);
}

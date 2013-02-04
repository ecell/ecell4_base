#define BOOST_TEST_MODULE "ODESimulator_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>

#include <ecell4/core/NetworkModel.hpp>
#include "../ODESimulator.hpp"

using namespace ecell4;
using namespace ecell4::ode;


BOOST_AUTO_TEST_CASE(ODESimulator_test_constructor)
{
    Real const volume(1e-18);

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    boost::shared_ptr<ODEWorld> world(new ODEWorld(volume));

    ODESimulator target(model, world);
}

BOOST_AUTO_TEST_CASE(ODESimulator_test_step)
{
    Real const volume(1e-18);

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    boost::shared_ptr<ODEWorld> world(new ODEWorld(volume));

    ODESimulator target(model, world);
    target.step();
}

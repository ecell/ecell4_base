#define BOOST_TEST_MODULE "SpatiocyteSimulator_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>

#include <ecell4/core/Position3.hpp>
#include <ecell4/core/NetworkModel.hpp>

#include "../SpatiocyteSimulator.hpp"

using namespace ecell4;
using namespace ecell4::spatiocyte;


BOOST_AUTO_TEST_CASE(SpatiocyteSimulator_test_constructor)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(1e-8);

    boost::shared_ptr<ecell4::Model> model(new ecell4::NetworkModel());
    boost::shared_ptr<SpatiocyteWorld> world(
        new SpatiocyteWorld(edge_lengths, voxel_radius));

    SpatiocyteSimulator target(model, world);
}

BOOST_AUTO_TEST_CASE(SpatiocyteSimulator_test_step)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(1e-8);

    boost::shared_ptr<ecell4::Model> model(new ecell4::NetworkModel());
    boost::shared_ptr<SpatiocyteWorld> world(
        new SpatiocyteWorld(edge_lengths, voxel_radius));

    SpatiocyteSimulator target(model, world);
    target.step();
}

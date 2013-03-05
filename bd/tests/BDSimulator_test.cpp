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

BOOST_AUTO_TEST_CASE(BDSimulator_test_step1)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    boost::shared_ptr<RandomNumberGenerator> rng(new GSLRandomNumberGenerator());

    boost::shared_ptr<Model> model(new NetworkModel());
    boost::shared_ptr<BDWorld> world(new BDWorld(edge_lengths, rng));

    BDSimulator target(model, world);
    target.step();
}

BOOST_AUTO_TEST_CASE(BDSimulator_test_step2)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    boost::shared_ptr<RandomNumberGenerator> rng(new GSLRandomNumberGenerator());

    boost::shared_ptr<Model> model(new NetworkModel());
    Species sp1("A", "2.5e-9", "1e-12");
    model->add_species(sp1);

    boost::shared_ptr<BDWorld> world(new BDWorld(edge_lengths, rng));
    world->new_particle(Particle(sp1, Position3(0, 0, 0), 2.5e-9, 1e-12));
    world->add_molecules(sp1, 10);

    BDSimulator target(model, world);
    target.step();
}

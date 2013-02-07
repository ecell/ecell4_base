#define BOOST_TEST_MODULE "GillespieSimulator_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/NetworkModel.hpp>

#include "../GillespieWorld.cpp"
#include "../GillespieSimulator.hpp"

using namespace ecell4;
using namespace ecell4::gillespie;

BOOST_AUTO_TEST_CASE(GillespieSimulator_test_step)
{
    Species sp1("A");
    Species sp2("B");

    ReactionRule rr1;
    rr1.set_k(5.0);
    rr1.add_reactant(sp1);
    rr1.add_product(sp2);

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_reaction_rule(rr1);
    Real vol(1.0);
    boost::shared_ptr<RandomNumberGenerator> rng(new GSLRandomNumberGenerator());
    rng->seed(time(NULL));
    boost::shared_ptr<GillespieWorld> world(new GillespieWorld(vol, rng));
    world->add_species(sp1);
    world->add_species(sp2);

    world->add_molecules(sp1, 10);
    world->add_molecules(sp2, 10);

    model->add_species(sp1);
    model->add_species(sp2);

    GillespieSimulator sim(model, world);

    sim.set_t(0.0);
    sim.step();

    BOOST_CHECK(0 < sim.t());
    BOOST_CHECK(world->num_molecules(sp1) == 9);

}


#define BOOST_TEST_MODULE "GillespieWorld_test"
//#define BOOST_TEST_NO_LIB

#include <boost/test/unit_test.hpp>

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/NetworkModel.hpp>

#include <ecell4/gillespie/GillespieWorld.cpp>
#include <ecell4/gillespie/GillespieSimulator.hpp>

using namespace ecell4;
using namespace ecell4::gillespie;

BOOST_AUTO_TEST_CASE(GillespieWorld_test)
{
    boost::shared_ptr<RandomNumberGenerator> rng(new GSLRandomNumberGenerator());
    rng->seed(time(NULL));
    Species sp1("A");
    Species sp2("B");

    const Real L(1.0);
    const Position3 edge_lengths(L, L, L);
    boost::shared_ptr<GillespieWorld> world(new GillespieWorld(edge_lengths, rng));

    world->add_molecules(sp1, 10);
    world->add_molecules(sp2, 20);
    world->set_t(0.5);

    BOOST_CHECK(world->t() == 0.5);
    BOOST_CHECK(world->num_molecules(sp1) == 10);
    BOOST_CHECK(world->num_molecules(sp2) == 20);
}

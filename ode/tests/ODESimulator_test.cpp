#define BOOST_TEST_MODULE "ODESimulator_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>

#include <ecell4/core/Species.hpp>
#include <ecell4/core/ReactionRule.hpp>
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

BOOST_AUTO_TEST_CASE(ODESimulator_test_step1)
{
    Real const volume(1e-18);

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    boost::shared_ptr<ODEWorld> world(new ODEWorld(volume));

    ODESimulator target(model, world);
    target.step(1.0);
}

BOOST_AUTO_TEST_CASE(ODESimulator_test_step2)
{
    Real const volume(1e-18);

    Species sp1("A"), sp2("B"), sp3("C");
    ReactionRule rr1;
    rr1.set_k(1.0);
    rr1.add_reactant(sp1);
    rr1.add_product(sp2);
    rr1.add_product(sp3);

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_species(sp1);
    model->add_species(sp2);
    model->add_species(sp3);
    model->add_reaction_rule(rr1);

    boost::shared_ptr<ODEWorld> world(new ODEWorld(volume));
    world->add_species(sp1);
    world->add_species(sp2);
    world->add_species(sp3);
    world->set_num_molecules(sp1, 60);

    ODESimulator target(model, world);

    // std::cout << target.t() << ":" << world->num_molecules(sp1)
    //           << ":" << world->num_molecules(sp2) << std::endl;
    target.step(1.0);
    // std::cout << target.t() << ":" << world->num_molecules(sp1)
    //           << ":" << world->num_molecules(sp2) << std::endl;

    // BOOST_ASSERT(false);
}

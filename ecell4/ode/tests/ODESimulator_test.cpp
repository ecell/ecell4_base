#define BOOST_TEST_MODULE "ODESimulator_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <ecell4/core/Species.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include "../ODESimulator.hpp"

using namespace ecell4;
using namespace ecell4::ode;


BOOST_AUTO_TEST_CASE(ODESimulator_test_constructor)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    boost::shared_ptr<ODEWorld> world(new ODEWorld(edge_lengths));

    ODESimulator target(model, world);
}

BOOST_AUTO_TEST_CASE(ODESimulator_test_step1)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    boost::shared_ptr<ODEWorld> world(new ODEWorld(edge_lengths));

    ODESimulator target(model, world);
    // target.step(1.0); //XXX: why not?
}

BOOST_AUTO_TEST_CASE(ODESimulator_test_step2)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);

    Species sp1("A"), sp2("B"), sp3("C");
    ReactionRule rr1;
    rr1.set_k(1.0);
    rr1.add_reactant(sp1);
    rr1.add_product(sp2);
    rr1.add_product(sp3);

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_species_attribute(sp1);
    model->add_species_attribute(sp2);
    model->add_species_attribute(sp3);
    model->add_reaction_rule(rr1);

    boost::shared_ptr<ODEWorld> world(new ODEWorld(edge_lengths));
    world->reserve_species(sp1);
    world->set_value(sp1, 60);

    ODESimulator target(model, world);

    // std::cout << target.t() << ":" << world->num_molecules(sp1)
    //           << ":" << world->num_molecules(sp2) << std::endl;
    target.step(1.0);
    // std::cout << target.t() << ":" << world->num_molecules(sp1)
    //           << ":" << world->num_molecules(sp2) << std::endl;

    // BOOST_ASSERT(false);
}

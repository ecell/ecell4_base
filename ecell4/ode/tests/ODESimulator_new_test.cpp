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
#include "../ODESimulator_New.hpp"


using namespace ecell4;
using namespace ecell4::ode;

BOOST_AUTO_TEST_CASE(ODESimulator_test_constructor)
{
    const Real L(1e-6);
    const Real3 edge_lengths(L, L, L);
    ODEWorld_New world(edge_lengths);
    BOOST_CHECK_EQUAL(world.t(), 0.);

    world.set_t(0.5);
    BOOST_CHECK_EQUAL(world.t(), 0.5);

    BOOST_CHECK_EQUAL(world.edge_lengths() , edge_lengths);
    BOOST_CHECK_EQUAL(world.volume(), edge_lengths[0] * edge_lengths[1] * edge_lengths[2]);

    // Compartment Space related
    Species sp1("A");
    Species sp2("B");
    Species sp3("C");
    world.add_molecules(sp1,  5.00);
    world.add_molecules(sp2, 10.00);
    world.add_molecules(sp3,  0.00);

    BOOST_CHECK_EQUAL(world.get_value(sp1), 5.0);
    BOOST_CHECK_EQUAL(world.get_value(sp2),10.0);
    BOOST_CHECK_EQUAL(world.get_value(sp3), 0.0);

    BOOST_CHECK_EQUAL(world.get_value_exact(sp1), 5.0);
    BOOST_CHECK_EQUAL(world.get_value_exact(sp2),10.0);
    BOOST_CHECK_EQUAL(world.get_value_exact(sp3), 0.0);

    BOOST_CHECK_EQUAL(world.has_species(sp1), true);
    BOOST_CHECK_EQUAL(world.has_species(Species("B")), true);
    BOOST_CHECK_EQUAL(world.has_species(Species("D")),false);
}

BOOST_AUTO_TEST_CASE(ODESimulator_test_ClasicRR)
{
    // Without ReactionDescriptor.
    const Real L(1e-6);
    const Real3 edge_lengths(L, L, L);
    std::shared_ptr<ODEWorld_New> world(new ODEWorld_New(edge_lengths));
    
    ReactionRule rr1;
    Species sp1("A"), sp2("B");
    rr1.add_reactant(sp1);
    rr1.add_product(sp2);
    rr1.set_k(1.0);

    Real sp1_initial_value = 100.;
    world->add_molecules(sp1, sp1_initial_value);

    std::shared_ptr<NetworkModel> new_model(new NetworkModel());
    {
        new_model->add_species_attribute(sp1);
        new_model->add_species_attribute(sp2);
        new_model->add_reaction_rule(rr1);
    }

    ODESolverType type = RUNGE_KUTTA_CASH_KARP54;
    ODESimulator_New target(world, new_model, type);

    Real dt = 0.1;
    target.set_dt(dt);

    BOOST_CHECK_EQUAL(target.t(), 0.);
    BOOST_CHECK_EQUAL(target.dt(), dt);

    BOOST_TEST(target.step(10) );   // must return true
    BOOST_TEST(world->get_value_exact(sp1) < sp1_initial_value);
}

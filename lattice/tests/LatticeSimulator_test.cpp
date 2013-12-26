#define BOOST_TEST_MODULE "LatticeSimulator_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "../LatticeSimulator.hpp"

using namespace ecell4;
using namespace ecell4::lattice;

BOOST_AUTO_TEST_CASE(LatticeSimulator_test_constructor)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(1e-8);
    const Integer N(60);

    const std::string D("1e-12"), radius("2.5e-9");

    ecell4::Species sp1("A", D, radius),
        sp2("B", D, radius),
        sp3("C", D, radius);
    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    (*model).add_species(sp1);
    (*model).add_species(sp2);
    (*model).add_species(sp3);

    boost::shared_ptr<GSLRandomNumberGenerator> rng;
    boost::shared_ptr<LatticeWorld> world(
            new LatticeWorld(edge_lengths, rng));

    LatticeSimulator sim(model, world);
}

BOOST_AUTO_TEST_CASE(LatticeSimulator_test_step_with_single_species)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(1e-8);
    const Integer N(60);

    const std::string D("1e-12"), radius("2.5e-9");

    ecell4::Species sp1("A", D, radius);
    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    (*model).add_species(sp1);

    boost::shared_ptr<GSLRandomNumberGenerator> rng;
    boost::shared_ptr<LatticeWorld> world(
            new LatticeWorld(edge_lengths, rng));

    //world->add_molecules(sp1, N / 2);
    world->add_molecule(sp1, 0);
    world->add_molecule(sp1, 2);

    //BOOST_ASSERT(world->num_molecules(sp1) == N / 2);
    BOOST_ASSERT(world->num_molecules(sp1) == 2);

    LatticeSimulator sim(model, world);

    //world->add_molecules(sp1, N / 2);
    //BOOST_ASSERT(world->num_molecules(sp1) == N);
    world->add_molecule(sp1, 4);
    BOOST_ASSERT(world->num_molecules(sp1) == 3);

    //sim.step();
}

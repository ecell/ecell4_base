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

    ecell4::Species sp1("A", radius, D),
        sp2("B", radius, D),
        sp3("C", radius, D);
    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    (*model).add_species(sp1);
    (*model).add_species(sp2);
    (*model).add_species(sp3);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
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

    ecell4::Species sp1("A", radius, D);
    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    (*model).add_species(sp1);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<LatticeWorld> world(
            new LatticeWorld(edge_lengths, rng));

    //world->add_molecules(sp1, N / 2);
    world->add_molecule(sp1, 161605);
    world->add_molecule(sp1, 161606);

    //BOOST_ASSERT(world->num_molecules(sp1) == N / 2);
    BOOST_ASSERT(world->num_molecules(sp1) == 2);

    LatticeSimulator sim(model, world);

    //world->add_molecules(sp1, N / 2);
    //BOOST_ASSERT(world->num_molecules(sp1) == N);
    world->add_molecule(sp1, 300000);
    BOOST_ASSERT(world->num_molecules(sp1) == 3);

    //sim.step();
}

BOOST_AUTO_TEST_CASE(LattiecSimulator_test_scheduler)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(1e-8);
    const Integer N(60);

    const std::string D1("1.0e-12"),
          D2("1.1e-12"),
          D3("1.2e-12"),
          radius("2.5e-9");

    ecell4::Species sp1("A", radius, D1),
        sp2("B", radius, D2),
        sp3("C", radius, D3);
    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    (*model).add_species(sp1);
    (*model).add_species(sp2);
    (*model).add_species(sp3);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<LatticeWorld> world(
            new LatticeWorld(edge_lengths, rng));

    Coord c1(658518), c2(300000), c3(486420);
    world->add_molecule(sp1, c1);
    world->add_molecule(sp2, c2);
    world->add_molecule(sp3, c3);

    LatticeSimulator sim(model, world);

    sim.initialize();

    sim.step();
    sim.step();
    sim.step();
    sim.step();
    sim.step();
    sim.step();
    sim.step();
    sim.step();
    sim.step();
    const MolecularTypeBase
        *mt1(world->get_molecular_type(sp1)),
        *mt2(world->get_molecular_type(sp2)),
        *mt3(world->get_molecular_type(sp3));
    std::vector<std::pair<Coord, ParticleID> >::const_iterator
        itr1(mt1->begin()),
        itr2(mt2->begin()),
        itr3(mt3->begin());
    c1 = (*itr1).first;
    c2 = (*itr2).first;
    c3 = (*itr3).first;

    sim.step();
    itr1 = mt1->begin();
    itr2 = mt2->begin();
    itr3 = mt3->begin();
    BOOST_ASSERT((*itr1).first == c1);
    BOOST_ASSERT((*itr2).first == c2);
    BOOST_ASSERT((*itr3).first != c3);
    c3 = (*itr3).first;

    sim.step();
    itr1 = mt1->begin();
    itr2 = mt2->begin();
    itr3 = mt3->begin();
    /*
    BOOST_ASSERT((*itr1).first == c1);
    BOOST_ASSERT((*itr2).first != c2);
    BOOST_ASSERT((*itr3).first == c3);
    */
    c2 = (*itr2).first;

    sim.step();
    itr1 = mt1->begin();
    itr2 = mt2->begin();
    itr3 = mt3->begin();
    /*
    BOOST_ASSERT((*itr1).first != c1);
    BOOST_ASSERT((*itr2).first == c2);
    BOOST_ASSERT((*itr3).first == c3);
    */
    c1 = (*itr1).first;

}


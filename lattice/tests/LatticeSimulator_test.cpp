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
    const Real voxel_radius(2.5e-9);

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
            new LatticeWorld(edge_lengths, voxel_radius, rng));

    LatticeSimulator sim(model, world);
}

BOOST_AUTO_TEST_CASE(LatticeSimulator_test_hdf5_save)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(2.5e-9);
    const Integer N(60);

    const std::string D("1e-12"), radius("2.5e-9");

    ecell4::Species sp("A", radius, D);
    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    (*model).add_species(sp);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<LatticeWorld> world(
            new LatticeWorld(edge_lengths, voxel_radius, rng));

    world->add_molecules(sp, N);
    BOOST_ASSERT(world->num_molecules(sp) == N);

    LatticeSimulator sim(model, world);

    H5::H5File fout("data.h5", H5F_ACC_TRUNC);
    const std::string hdf5path("/");
    world->save(&fout, hdf5path);
}

BOOST_AUTO_TEST_CASE(LatticeSimulator_test_step_with_single_particle)
{
    const Real L(2.5e-8);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(2.5e-9);

    const std::string D("1e-12"), radius("2.5e-9");

    ecell4::Species sp("A", radius, D);
    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    (*model).add_species(sp);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<LatticeWorld> world(
            new LatticeWorld(edge_lengths, voxel_radius, rng));

    BOOST_CHECK(world->add_molecule(sp, 36));

    LatticeSimulator sim(model, world);

    const std::string hdf5path("/");

    for (int i(0); i < 50; ++i)
    {
        std::ostringstream oss;
        oss << "data_with_single_particle_";
        if (i < 10)
        {
            oss << "0" << i;
        }
        else
        {
            oss << i;
        }
        oss << ".h5";
        H5::H5File fout(oss.str(), H5F_ACC_TRUNC);
        sim.step();
        world->save(&fout, hdf5path);
    }
}

BOOST_AUTO_TEST_CASE(LatticeSimulator_test_step_with_single_species)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(2.5e-9);
    const Integer N(60);

    const std::string D("1e-12"), radius("2.5e-9");

    ecell4::Species sp("A", radius, D);
    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    (*model).add_species(sp);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<LatticeWorld> world(
            new LatticeWorld(edge_lengths, voxel_radius, rng));

    world->add_molecules(sp, N / 2);

    BOOST_ASSERT(world->num_molecules(sp) == N / 2);

    LatticeSimulator sim(model, world);

    world->add_molecules(sp, N / 2);
    BOOST_ASSERT(world->num_molecules(sp) == N);

    sim.step();
}

BOOST_AUTO_TEST_CASE(LatticeSimulator_test_save_step_with_single_species)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(2.5e-9);
    const Integer N(60);

    const std::string D("1e-12"), radius("2.5e-9");

    ecell4::Species sp("A", radius, D);
    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    (*model).add_species(sp);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<LatticeWorld> world(
            new LatticeWorld(edge_lengths, voxel_radius, rng));

    LatticeSimulator sim(model, world);

    world->add_molecules(sp, N);

    const std::string hdf5path("/");

    for (int i(0); i < 50; ++i)
    {
        std::ostringstream oss;
        oss << "data_with_single_species_";
        if (i < 10)
        {
            oss << "0" << i;
        }
        else
        {
            oss << i;
        }
        oss << ".h5";
        H5::H5File fout(oss.str(), H5F_ACC_TRUNC);
        sim.step();
        world->save(&fout, hdf5path);
    }
}

BOOST_AUTO_TEST_CASE(LatticeSimulator_test_save_step_with_periodic)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(2.5e-9);
    const Integer N(60);

    const std::string D("1e-12"), radius("2.5e-9");

    ecell4::Species sp("A", radius, D);
    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    (*model).add_species(sp);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<LatticeWorld> world(
            new LatticeWorld(edge_lengths, voxel_radius, rng));

    LatticeSimulator sim(model, world);

    world->add_molecules(sp, N);

    const std::string hdf5path("/");

    for (int i(0); i < 50; ++i)
    {
        std::ostringstream oss;
        oss << "data_with_single_species_";
        if (i < 10)
        {
            oss << "0" << i;
        }
        else
        {
            oss << i;
        }
        oss << ".h5";
        H5::H5File fout(oss.str(), H5F_ACC_TRUNC);
        sim.step();
        world->save(&fout, hdf5path);
    }
}

BOOST_AUTO_TEST_CASE(LatticeSimulator_test_unimolecular_reaction)
{
    const Real L(2.5e-8);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(2.5e-9);
    const std::string radius("1.25e-9");
    const ecell4::Species sp1("A", radius, "1.0e-12"),
          sp2("B", radius, "1.1e-12"),
          sp3("C", "2.5e-9", "1.2e-12");

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_species(sp1);
    model->add_species(sp2);
    model->add_species(sp3);

    model->add_reaction_rule(create_unimolecular_reaction_rule(sp1,sp3,1e6));

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<LatticeWorld> world(
            new LatticeWorld(edge_lengths, voxel_radius, rng));

    LatticeSimulator sim(model, world);

    BOOST_CHECK(world->add_molecules(sp1, 25));
    BOOST_CHECK(world->add_molecules(sp2, 25));

    H5::H5File fout_before("data_unimolecular_reaction_single0.h5", H5F_ACC_TRUNC);
    world->save(&fout_before, "/");
    for (Integer i(0); i < 10; ++i)
    {
        sim.step();
    }
    BOOST_ASSERT(world->num_molecules(sp3) > 0);
    BOOST_ASSERT(25 - world->num_molecules(sp1) == world->num_molecules(sp3));
    H5::H5File fout_after("data_unimolecular_reaction_single1.h5", H5F_ACC_TRUNC);
    world->save(&fout_after, "/");
}

BOOST_AUTO_TEST_CASE(LatticeSimulator_test_binding_reaction)
{
    const Real L(2.5e-8);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(2.5e-9);
    const std::string radius("1.25e-9");
    const ecell4::Species sp1("A", radius, "1.0e-12"),
          sp2("B", radius, "1.1e-12"),
          sp3("C", "2.5e-9", "1.2e-12");

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_species(sp1);
    model->add_species(sp2);
    model->add_species(sp3);

    model->add_reaction_rule(create_binding_reaction_rule(sp1,sp2,sp3,5.));

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<LatticeWorld> world(
            new LatticeWorld(edge_lengths, voxel_radius, rng));

    LatticeSimulator sim(model, world);

    BOOST_CHECK(world->add_molecules(sp1, 25));
    BOOST_CHECK(world->add_molecules(sp2, 25));

    H5::H5File fout_before("data_binging_reaction0.h5", H5F_ACC_TRUNC);
    world->save(&fout_before, "/");
    sim.step();
    sim.step();
    BOOST_ASSERT(world->num_molecules(sp3) > 0);
    BOOST_ASSERT(25 - world->num_molecules(sp1) == world->num_molecules(sp3));
    BOOST_ASSERT(25 - world->num_molecules(sp2) == world->num_molecules(sp3));
    H5::H5File fout_after("data_binding_reaction1.h5", H5F_ACC_TRUNC);
    world->save(&fout_after, "/");
}

BOOST_AUTO_TEST_CASE(LatticeSimulator_test_unbinding_reaction)
{
    const Real L(2.5e-8);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(2.5e-9);
    const std::string radius("1.25e-9");
    const ecell4::Species sp1("A", radius, "1.0e-12"),
          sp2("B", radius, "1.1e-12"),
          sp3("C", "2.5e-9", "1.2e-12");

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_species(sp1);
    model->add_species(sp2);
    model->add_species(sp3);

    model->add_reaction_rule(create_unbinding_reaction_rule(sp1,sp2,sp3,1e5));

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<LatticeWorld> world(
            new LatticeWorld(edge_lengths, voxel_radius, rng));

    LatticeSimulator sim(model, world);

    BOOST_CHECK(world->add_molecules(sp1, 25));

    H5::H5File fout_before("data_unbinding_reaction0.h5", H5F_ACC_TRUNC);
    world->save(&fout_before, "/");
    for (Integer i(0); i < 10; ++i)
    {
        sim.step();
    }
    BOOST_ASSERT(world->num_molecules(sp1) < 25);
    BOOST_ASSERT(25 - world->num_molecules(sp1) == world->num_molecules(sp2));
    BOOST_ASSERT(25 - world->num_molecules(sp1) == world->num_molecules(sp3));
    H5::H5File fout_after("data_unbinding_reaction1.h5", H5F_ACC_TRUNC);
    world->save(&fout_after, "/");
}

BOOST_AUTO_TEST_CASE(LatticeSimulator_test_degradation_reaction)
{
    const Real L(2.5e-8);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(2.5e-9);
    const std::string radius("1.25e-9");
    const ecell4::Species sp1("A", radius, "1.0e-12");

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_species(sp1);

    model->add_reaction_rule(create_degradation_reaction_rule(sp1,1e5));

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<LatticeWorld> world(
            new LatticeWorld(edge_lengths, voxel_radius, rng));

    LatticeSimulator sim(model, world);

    BOOST_CHECK(world->add_molecules(sp1, 25));

    H5::H5File fout_before("data_degradation_reaction0.h5", H5F_ACC_TRUNC);
    world->save(&fout_before, "/");
    for (Integer i(0); i < 10; ++i)
    {
        sim.step();
    }
    BOOST_ASSERT(world->num_molecules(sp1) < 25);
    H5::H5File fout_after("data_degradation_reaction1.h5", H5F_ACC_TRUNC);
    world->save(&fout_after, "/");
}

BOOST_AUTO_TEST_CASE(LattiecSimulator_test_scheduler)
{
    std::cout << " <<LattiecSimulator_test_scheduler>> ";
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(2.5e-9);

    const std::string D1("1.0e-12"),
          D2("1.1e-12"),
          D3("1.2e-12"),
          radius("2.5e-9");

    const ecell4::Species sp1("A", radius, D1),
        sp2("B", radius, D2),
        sp3("C", radius, D3);
    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    (*model).add_species(sp1);
    (*model).add_species(sp2);
    (*model).add_species(sp3);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<LatticeWorld> world(
            new LatticeWorld(edge_lengths, voxel_radius, rng));

    Coord c1(world->global2coord(Global(40,34,56))),
          c2(world->global2coord(Global(32,50,24))),
          c3(world->global2coord(Global(60,36,89)));
    BOOST_CHECK(world->add_molecule(sp1, c1));
    BOOST_CHECK(world->add_molecule(sp2, c2));
    BOOST_CHECK(world->add_molecule(sp3, c3));

    LatticeSimulator sim(model, world);

    sim.initialize();

    const MolecularTypeBase
        *mt1(world->get_molecular_type(sp1)),
        *mt2(world->get_molecular_type(sp2)),
        *mt3(world->get_molecular_type(sp3));
    std::vector<std::pair<Coord, ParticleID> >::const_iterator
        itr1(mt1->begin()),
        itr2(mt2->begin()),
        itr3(mt3->begin());

    BOOST_ASSERT(itr1 != mt1->end());
    BOOST_ASSERT(itr2 != mt2->end());
    BOOST_ASSERT(itr3 != mt3->end());

    c1 = (*itr1).first;
    c2 = (*itr2).first;
    c3 = (*itr3).first;

    sim.step();
    itr1 = mt1->begin();
    itr2 = mt2->begin();
    itr3 = mt3->begin();
    std::cout << "<itr1: " << (*itr1).first << "> ";
    std::cout << "<itr2: " << (*itr2).first << "> ";
    std::cout << "<itr3: " << (*itr3).first << "> ";
    BOOST_ASSERT((*itr1).first == c1);
    BOOST_ASSERT((*itr2).first == c2);
    BOOST_ASSERT((*itr3).first != c3);
    c3 = (*itr3).first;

    sim.step();
    itr1 = mt1->begin();
    itr2 = mt2->begin();
    itr3 = mt3->begin();
    BOOST_ASSERT((*itr1).first == c1);
    BOOST_ASSERT((*itr2).first != c2);
    BOOST_ASSERT((*itr3).first == c3);
    c2 = (*itr2).first;

    sim.step();
    itr1 = mt1->begin();
    itr2 = mt2->begin();
    itr3 = mt3->begin();
    BOOST_ASSERT((*itr1).first != c1);
    BOOST_ASSERT((*itr2).first == c2);
    BOOST_ASSERT((*itr3).first == c3);
    c1 = (*itr1).first;

}


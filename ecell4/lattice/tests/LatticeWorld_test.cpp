#define BOOST_TEST_MODULE "LatticeWorld_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/floating_point_comparison.hpp>

#include "../LatticeWorld.hpp"

using namespace ecell4;
using namespace ecell4::lattice;

const Real DEFAULT_VOXEL_RADIUS = 1e-8;

BOOST_AUTO_TEST_CASE(LatticeWorld_test_constructor)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    const Position3 edge_lengths(1e-6, 1e-6, 1e-6);
    const Real voxel_radius(DEFAULT_VOXEL_RADIUS);
    LatticeWorld world(edge_lengths, voxel_radius, rng);
}

BOOST_AUTO_TEST_CASE(LatticeWorld_test_t)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    const Position3 edge_lengths(1e-6, 1e-6, 1e-6);
    const Real voxel_radius(DEFAULT_VOXEL_RADIUS);
    LatticeWorld world(edge_lengths, voxel_radius, rng);
    BOOST_CHECK_EQUAL(world.t(), 0);
    world.set_t(23.4);
    BOOST_CHECK_EQUAL(world.t(), 23.4);
}

BOOST_AUTO_TEST_CASE(LatticeWorld_test_num_species)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    const Position3 edge_lengths(1e-6, 1e-6, 1e-6);
    const Real voxel_radius(DEFAULT_VOXEL_RADIUS);
    LatticeWorld world(edge_lengths, voxel_radius, rng);
    BOOST_CHECK_EQUAL(world.num_species(), 0);
}

BOOST_AUTO_TEST_CASE(LatticeWorld_test_has_species)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    const Position3 edge_lengths(1e-6, 1e-6, 1e-6);
    const Real voxel_radius(DEFAULT_VOXEL_RADIUS);
    LatticeWorld world(edge_lengths, voxel_radius, rng);
    Species sp(std::string("Species"));
    BOOST_CHECK(!world.has_species(sp));
}

BOOST_AUTO_TEST_CASE(LatticeWorld_test_list_particles)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    const Position3 edge_lengths(1e-6, 1e-6, 1e-6);
    const Real voxel_radius(DEFAULT_VOXEL_RADIUS);
    LatticeWorld world(edge_lengths, voxel_radius, rng);
    std::vector<std::pair<ParticleID, Particle> > particles(world.list_particles());
    BOOST_CHECK_EQUAL(particles.size(), 0);
}

BOOST_AUTO_TEST_CASE(LatticeWorld_test_update_particles)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    SerialIDGenerator<ParticleID> sidgen;
    const Position3 edge_lengths(1e-6, 1e-6, 1e-6);
    const Real voxel_radius(DEFAULT_VOXEL_RADIUS);
    LatticeWorld world(edge_lengths, voxel_radius, rng);

    ParticleID pid(sidgen());
    Species sp(std::string("A"));
    const Position3 pos(2e-7, 1e-7, 0);
    Real r(0);
    Real d(0);
    Particle p(sp, pos, r, d);

    world.update_particle(pid, p);

    BOOST_CHECK(world.has_species(sp));
    BOOST_CHECK(world.has_particle(pid));
    BOOST_CHECK_EQUAL(world.list_particles().size(), 1);
    BOOST_CHECK_EQUAL(world.list_particles(sp).size(), 1);
}

// BOOST_AUTO_TEST_CASE(LatticeWorld_test_register_species)
// {
//     const Position3 edge_lengths(1e-6,1e-6,1e-6);
//     const Real voxel_radius(DEFAULT_VOXEL_RADIUS);
//     boost::shared_ptr<GSLRandomNumberGenerator>
//         rng(new GSLRandomNumberGenerator());
//     LatticeWorld world(edge_lengths, voxel_radius, rng);
// 
//     Species sp(std::string("TEST"));
// 
//     BOOST_CHECK(world.register_species(sp));
//     BOOST_CHECK(world.has_species(sp));
// 
//     std::vector<Species> list;
//     list.push_back(sp);
// 
//     BOOST_CHECK(list == world.list_species());
// }

BOOST_AUTO_TEST_CASE(LatticeWorld_test_add_molecule)
{
    const Position3 edge_lengths(1e-6,1e-6,1e-6);
    const Real voxel_radius(2.5e-9);
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    LatticeWorld world(edge_lengths, voxel_radius, rng);

    Species sp(std::string("TEST"));
    sp.set_attribute("radius", "2.5e-9");
    sp.set_attribute("D", "1e-12");

    LatticeWorld::private_coordinate_type coord(486420);
    BOOST_CHECK(world.place_voxel_private(sp, coord).second);
    BOOST_CHECK_EQUAL(world.num_particles(sp), 1);

    MolecularTypeBase* mt(world.get_molecular_type_private(coord));
    BOOST_CHECK(!mt->is_vacant());
}

BOOST_AUTO_TEST_CASE(LatticeWorld_test_add_molecules)
{
    const Position3 edge_lengths(1e-6,1e-6,1e-6);
    const Real voxel_radius(DEFAULT_VOXEL_RADIUS);
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    LatticeWorld world(edge_lengths, voxel_radius, rng);

    Species sp(std::string("TEST"));
    sp.set_attribute("radius", "2.5e-9");
    sp.set_attribute("D", "1e-12");
    const Integer N(60);

    BOOST_CHECK(world.add_molecules(sp, N));
    BOOST_CHECK_EQUAL(world.num_particles(sp), N);
}

BOOST_AUTO_TEST_CASE(LatticeWorld_test_move)
{
    const Position3 edge_lengths(1e-6,1e-6,1e-6);
    const Real voxel_radius(2.5e-9);
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    LatticeWorld world(edge_lengths, voxel_radius, rng);

    Species sp(std::string("TEST"));
    sp.set_attribute("radius", "2.5e-9");
    sp.set_attribute("D", "1e-12");

    LatticeWorld::coordinate_type from(1034), to(786420);

    LatticeWorld::private_coordinate_type private_from(
            world.coord2private(from));
    BOOST_CHECK(world.place_voxel_private(sp, private_from).second);

    LatticeWorld::private_coordinate_type private_to(
            world.coord2private(to));
    BOOST_CHECK(world.move(from, to));

    MolecularTypeBase* mt(world.get_molecular_type_private(private_to));
    BOOST_CHECK(!mt->is_vacant());

    BOOST_CHECK(world.move(from, to));
}

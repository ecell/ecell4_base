#define BOOST_TEST_MODULE "LatticeWorld_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "../LatticeWorld.hpp"
#include <ecell4/core/SerialIDGenerator.hpp>

using namespace ecell4;
using namespace ecell4::lattice;


BOOST_AUTO_TEST_CASE(LatticeWorld_test_constructor)
{
    boost::shared_ptr<GSLRandomNumberGenerator> rng;
    Position3 edge_lengths(1e-6, 1e-6, 1e-6);
    LatticeWorld world(edge_lengths, rng);
}

BOOST_AUTO_TEST_CASE(LatticeWorld_test_t)
{
    boost::shared_ptr<GSLRandomNumberGenerator> rng;
    Position3 edge_lengths(1e-6, 1e-6, 1e-6);
    LatticeWorld world(edge_lengths, rng);
    BOOST_CHECK_EQUAL(world.t(), 0);
    world.set_t(23.4);
    BOOST_CHECK_EQUAL(world.t(), 23.4);
}

BOOST_AUTO_TEST_CASE(LatticeWorld_test_num_species)
{
    boost::shared_ptr<GSLRandomNumberGenerator> rng;
    Position3 edge_lengths(1e-6, 1e-6, 1e-6);
    LatticeWorld world(edge_lengths, rng);
    BOOST_CHECK_EQUAL(world.num_species(), 0);
}

BOOST_AUTO_TEST_CASE(LatticeWorld_test_has_species)
{
    boost::shared_ptr<GSLRandomNumberGenerator> rng;
    Position3 edge_lengths(1e-6, 1e-6, 1e-6);
    LatticeWorld world(edge_lengths, rng);
    Species sp(std::string("Species"));
    BOOST_CHECK(!world.has_species(sp));
}

BOOST_AUTO_TEST_CASE(LatticeWorld_test_list_particles)
{
    boost::shared_ptr<GSLRandomNumberGenerator> rng;
    Position3 edge_lengths(1e-6, 1e-6, 1e-6);
    LatticeWorld world(edge_lengths, rng);
    std::vector<std::pair<ParticleID, Particle> > particles(world.list_particles());
    BOOST_CHECK_EQUAL(particles.size(), 0);
}

BOOST_AUTO_TEST_CASE(LatticeWorld_test_update_particles)
{
    boost::shared_ptr<GSLRandomNumberGenerator> rng;
    Position3 edge_lengths(1e-6, 1e-6, 1e-6);
    LatticeWorld world(edge_lengths, rng);
    SerialIDGenerator<ParticleID> sidgen;

    ParticleID pid(sidgen());
    Species sp(std::string("Species A"));
    Position3 pos(2e-7, 1e-7, 0);
    Real r(0);
    Real d(0);
    Particle p(sp, pos, r, d);

    world.update_particle(pid, p);

    BOOST_CHECK(world.has_species(sp));
    BOOST_CHECK(world.has_particle(pid));
    BOOST_CHECK_EQUAL(world.list_particles().size(), 1);
    BOOST_CHECK_EQUAL(world.list_particles(sp).size(), 1);
}

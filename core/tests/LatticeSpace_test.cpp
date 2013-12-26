#define BOOST_TEST_MODULE "ReactionRule_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "../LatticeSpace.hpp"
#include "../SerialIDGenerator.hpp"

using namespace ecell4;

BOOST_AUTO_TEST_CASE(LatticeSpace_test_constructor)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_num_species)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);
    BOOST_CHECK_EQUAL(lspace.num_species(), 0);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_has_species)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);
    const Species &sp = Species("TEST");
    BOOST_CHECK(!lspace.has_species(sp));
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_update_particle)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);

    SerialIDGenerator<ParticleID> sidgen;
    ParticleID id(sidgen());

    Species sp(std::string("TEST"));
    Position3 pos(2e-7, 0, 0);
    Real r(1.0);
    Real d(2.3);
    Particle particle(sp, pos, r, d);

    BOOST_CHECK(lspace.update_particle(id, particle));
    BOOST_CHECK(lspace.has_species(sp));
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_num_particles)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);
    BOOST_CHECK_EQUAL(lspace.num_particles(), 0);

    SerialIDGenerator<ParticleID> sidgen;
    ParticleID id(sidgen());
    Species sp(std::string("TEST"));
    Position3 pos(2e-7, 0, 0);
    Real r(1.0);
    Real d(2.3);
    Particle particle(sp, pos, r, d);

    ParticleID a_id(sidgen());
    Species a(std::string("ANOTHER"));
    Position3 pos1(1e-7, 2e-7, 0);
    Real r1(1.1);
    Real d1(4.3);
    Particle another(a, pos1, r1, d1);

    BOOST_CHECK(lspace.update_particle(id, particle));
    BOOST_CHECK(lspace.update_particle(a_id, another));
    BOOST_CHECK_EQUAL(lspace.num_particles(sp), 1);
    BOOST_CHECK_EQUAL(lspace.num_particles(), 2);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_list_particles)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);

    SerialIDGenerator<ParticleID> sidgen;
    ParticleID id(sidgen());
    Species sp(std::string("TEST"));
    Position3 pos(2e-7, 0, 0);
    Real r(1.0);
    Real d(2.3);
    Particle particle(sp, pos, r, d);

    ParticleID a_id(sidgen());
    Species a(std::string("ANOTHER"));
    Position3 pos1(1e-7, 2e-7, 0);
    Real r1(1.1);
    Real d1(4.3);
    Particle another(a, pos1, r1, d1);

    BOOST_CHECK(lspace.update_particle(id, particle));
    BOOST_CHECK(lspace.update_particle(a_id, another));

    typedef std::vector<std::pair<ParticleID, Particle> > vector;

    vector test_list(lspace.list_particles(sp));
    vector list(lspace.list_particles());
    BOOST_CHECK_EQUAL(list.size(), 2);
    BOOST_CHECK_EQUAL(test_list.size(), 1);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_add_species)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);

    Species sp(std::string("TEST"));

    BOOST_CHECK(lspace.add(sp));
    BOOST_CHECK(lspace.has_species(sp));

    std::vector<Species> list;
    list.push_back(sp);

    BOOST_CHECK(list == lspace.list_species());
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_add)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);
    SerialIDGenerator<ParticleID> sidgen;

    Species sp(std::string("TEST"));
    BOOST_CHECK(lspace.add(sp));

    Coord coord(2);
    ParticleID pid(sidgen());
    BOOST_CHECK(lspace.add(sp, coord, pid));
    BOOST_CHECK_EQUAL(lspace.num_particles(sp), 1);

    MolecularTypeBase* mt(lspace.get_molecular_type(coord));
    BOOST_CHECK(!mt->is_vacant());
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_move)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);
    SerialIDGenerator<ParticleID> sidgen;

    Species sp(std::string("TEST"));
    BOOST_CHECK(lspace.add(sp));

    Coord coord(2);
    ParticleID pid(sidgen());
    BOOST_CHECK(lspace.add(sp, coord, pid));

    Coord to_coord(1);
    BOOST_CHECK(lspace.move(coord, to_coord));

    MolecularTypeBase* mt(lspace.get_molecular_type(to_coord));
    BOOST_CHECK(!mt->is_vacant());

    BOOST_CHECK(!lspace.move(coord, to_coord));
}


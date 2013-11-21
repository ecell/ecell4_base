#define BOOST_TEST_MODULE "ReactionRule_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "../LatticeSpace.hpp"
#include "../SerialIDGenerator.hpp"

using namespace ecell4;


BOOST_AUTO_TEST_CASE(LatticeSpace_test_constructor)
{
    LatticeSpace lspace;
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_num_species)
{
    LatticeSpace lspace;
    BOOST_CHECK_EQUAL(lspace.num_species(), 0);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_has_species)
{
    LatticeSpace lspace;
    const Species &sp = Species("TEST");
    BOOST_CHECK(!lspace.has_species(sp));
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_update_sparticle)
{
    LatticeSpace lspace;
    SerialIDGenerator<ParticleID> sidgen;
    ParticleID id(sidgen());
    Species sp(std::string("TEST"));
    Position3 pos(2, 0, 0);
    Real r(1.0);
    Real d(2.3);
    Particle particle(sp, pos, r, d);
    BOOST_CHECK(lspace.update_particle(id, particle));
    BOOST_CHECK(lspace.has_species(sp));
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_num_particles)
{
    LatticeSpace lspace;
    BOOST_CHECK_EQUAL(lspace.num_particles(), 0);

    SerialIDGenerator<ParticleID> sidgen;
    ParticleID id(sidgen());
    Species sp(std::string("TEST"));
    Position3 pos(2, 0, 0);
    Real r(1.0);
    Real d(2.3);
    Particle particle(sp, pos, r, d);

    ParticleID a_id(sidgen());
    Species a(std::string("ANOTHER"));
    Position3 pos1(1, 2, 0);
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
    LatticeSpace lspace;

    SerialIDGenerator<ParticleID> sidgen;
    ParticleID id(sidgen());
    Species sp(std::string("TEST"));
    Position3 pos(2, 0, 0);
    Real r(1.0);
    Real d(2.3);
    Particle particle(sp, pos, r, d);

    ParticleID a_id(sidgen());
    Species a(std::string("ANOTHER"));
    Position3 pos1(1, 2, 0);
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

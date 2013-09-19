#define BOOST_TEST_MODULE "ReactionRule_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "../LatticeSpace.hpp"

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
    const Species sp = Species("TEST");
    BOOST_CHECK(!lspace.has_species(sp));
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_num_particles)
{
    LatticeSpace lspace;
    BOOST_CHECK_EQUAL(lspace.num_particles(), 0);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_update_sparticle)
{
    LatticeSpace lspace;
    ParticleID id = ParticleID();
    Species sp = Species("TEST");
    SParticle sparticle = {
        .coord = 6,
        .species = sp,
    };
    BOOST_CHECK(lspace.update_sparticle(id, sparticle));
    BOOST_CHECK(lspace.has_species(sp));
}

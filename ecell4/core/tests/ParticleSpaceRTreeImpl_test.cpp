#define BOOST_TEST_MODULE "ParticleSpaceRTreeImpl_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/tools/floating_point_comparison.hpp>

#include <ecell4/core/ParticleSpaceRTreeImpl.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>

using namespace ecell4;

// a set of parameters used in the test cases below
struct Fixture
{
    typedef ParticleSpaceRTreeImpl particle_space_type;

    const Real3 edge_lengths;
    const Real radius;

    Fixture() :
        edge_lengths(1, 1, 1),
        radius(0.005)
    {}
};

BOOST_FIXTURE_TEST_SUITE(suite, Fixture) // submit parameters

BOOST_AUTO_TEST_CASE(ParticleSpace_test_constructor)
{
    std::unique_ptr<ParticleSpace> space(new particle_space_type(edge_lengths));

    BOOST_CHECK_EQUAL((*space).list_species().size(), 0);
    BOOST_CHECK_EQUAL((*space).num_particles(), 0);
    BOOST_CHECK_EQUAL((*space).edge_lengths(), edge_lengths);
    BOOST_CHECK_EQUAL((*space).volume(), 1.0);
    BOOST_CHECK_EQUAL((*space).t(), 0.0);
}

BOOST_AUTO_TEST_CASE(ParticleSpace_test_update_remove)
{
    std::unique_ptr<ParticleSpace> space(new particle_space_type(edge_lengths));
    SerialIDGenerator<ParticleID> pidgen;

    const ParticleID pid1 = pidgen();
    const Species sp1 = Species("A");
    const Species sp2 = Species("B");

    BOOST_CHECK((*space).update_particle(pid1, Particle(sp1, edge_lengths * 0.5, radius, 0)));
    BOOST_CHECK_EQUAL((*space).num_particles(), 1);
    BOOST_CHECK_EQUAL((*space).num_particles(sp1), 1);
    BOOST_CHECK_EQUAL((*space).num_particles(sp2), 0);

    BOOST_CHECK(! (*space).update_particle(pid1, Particle(sp1, edge_lengths * 0.25, radius, 0)));
    BOOST_CHECK_EQUAL((*space).num_particles(), 1);
    BOOST_CHECK_EQUAL((*space).num_particles(sp1), 1);
    BOOST_CHECK_EQUAL((*space).num_particles(sp2), 0);

    BOOST_CHECK(! (*space).update_particle(pid1, Particle(sp2, edge_lengths * 0.1, radius, 0)));
    BOOST_CHECK_EQUAL((*space).num_particles(), 1);
    BOOST_CHECK_EQUAL((*space).num_particles(sp1), 0);
    BOOST_CHECK_EQUAL((*space).num_particles(sp2), 1);

    (*space).remove_particle(pid1);
    BOOST_CHECK_EQUAL((*space).num_particles(), 0);
    BOOST_CHECK_EQUAL((*space).num_particles(sp1), 0);
    BOOST_CHECK_EQUAL((*space).num_particles(sp2), 0);
}

BOOST_AUTO_TEST_CASE(ParticleSpace_test_remove)
{
    std::unique_ptr<ParticleSpace> space(new particle_space_type(edge_lengths));
    SerialIDGenerator<ParticleID> pidgen;

    const ParticleID pid1 = pidgen();
    const ParticleID pid2 = pidgen();
    const ParticleID pid3 = pidgen();
    const Species sp1 = Species("A");

    BOOST_CHECK((*space).update_particle(pid1, Particle(sp1, edge_lengths * 0.5, radius, 0)));
    BOOST_CHECK((*space).update_particle(pid2, Particle(sp1, edge_lengths * 0.25, radius, 0)));
    BOOST_CHECK((*space).update_particle(pid3, Particle(sp1, edge_lengths * 0.75, radius, 0)));
    BOOST_CHECK_EQUAL((*space).num_particles(sp1), 3);

    (*space).remove_particle(pid2);
    BOOST_CHECK_EQUAL((*space).num_particles(sp1), 2);
    (*space).remove_particle(pid3);
    BOOST_CHECK_EQUAL((*space).num_particles(sp1), 1);
    (*space).remove_particle(pid1);
    BOOST_CHECK_EQUAL((*space).num_particles(sp1), 0);
}

BOOST_AUTO_TEST_CASE(ParticleSpace_test_exception)
{
    std::unique_ptr<ParticleSpace> space(new particle_space_type(edge_lengths));
    SerialIDGenerator<ParticleID> pidgen;
    const ParticleID pid1 = pidgen();

    BOOST_CHECK_THROW((*space).remove_particle(pid1), NotFound);
}

BOOST_AUTO_TEST_CASE(ParticleSpace_test_list_particles_within_radius)
{
    std::unique_ptr<ParticleSpace> space(new particle_space_type(edge_lengths));
    SerialIDGenerator<ParticleID> pidgen;

    const ParticleID pid1 = pidgen();
    const ParticleID pid2 = pidgen();
    const Species sp1 = Species("A");

    BOOST_CHECK((*space).update_particle(pid1, Particle(sp1, Real3(0.5, 0.5, 0.5), radius, 0)));
    BOOST_CHECK_EQUAL((*space).num_particles(sp1), 1);

    {
        const std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
            retval = (*space).list_particles_within_radius(Real3(0.509, 0.5, 0.5), radius);
        BOOST_CHECK_EQUAL(retval.size(), 1);
        BOOST_CHECK_EQUAL(retval[0].first.first, pid1);
        BOOST_CHECK_CLOSE(retval[0].second, 0.009 - radius, 1e-6);
    }
    BOOST_CHECK_EQUAL((*space).list_particles_within_radius(Real3(0.509, 0.5, 0.5), radius, pid1).size(), 0);
    BOOST_CHECK_EQUAL((*space).list_particles_within_radius(Real3(0.511, 0.5, 0.5), radius).size(), 0);

    BOOST_CHECK((*space).update_particle(pid2, Particle(sp1, Real3(0.511, 0.5, 0.5), radius, 0)));
    BOOST_CHECK_EQUAL((*space).num_particles(sp1), 2);

    BOOST_CHECK_EQUAL((*space).list_particles_within_radius(Real3(0.509, 0.5, 0.5), radius).size(), 2);
    (*space).remove_particle(pid1);
    {
        const std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
            retval = (*space).list_particles_within_radius(Real3(0.509, 0.5, 0.5), radius);
        BOOST_CHECK_EQUAL(retval.size(), 1);
        BOOST_CHECK_EQUAL(retval[0].first.first, pid2);
        BOOST_CHECK_CLOSE(retval[0].second, 0.002 - radius, 1e-6);
    }
}

BOOST_AUTO_TEST_CASE(ParticleSpace_test_list_particles_within_radius_boundary_x)
{
    std::unique_ptr<ParticleSpace> space(new particle_space_type(edge_lengths));
    SerialIDGenerator<ParticleID> pidgen;

    const ParticleID pid1 = pidgen();
    const ParticleID pid2 = pidgen();
    const Species sp1 = Species("A");

    BOOST_CHECK((*space).update_particle(pid1, Particle(sp1, Real3(0.005, 0.5, 0.5), radius, 0)));
    BOOST_CHECK_EQUAL((*space).num_particles(sp1), 1);

    {
        const std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
            retval = (*space).list_particles_within_radius(Real3(0.996, 0.5, 0.5), radius);
        BOOST_CHECK_EQUAL(retval.size(), 1);
        BOOST_CHECK_EQUAL(retval[0].first.first, pid1);
        BOOST_CHECK_CLOSE(retval[0].second, 0.009 - radius, 1e-6);
    }
    BOOST_CHECK_EQUAL((*space).list_particles_within_radius(Real3(0.996, 0.5, 0.5), radius, pid1).size(), 0);
    BOOST_CHECK_EQUAL((*space).list_particles_within_radius(Real3(0.994, 0.5, 0.5), radius).size(), 0);

    BOOST_CHECK((*space).update_particle(pid2, Particle(sp1, Real3(0.994, 0.5, 0.5), radius, 0)));
    BOOST_CHECK_EQUAL((*space).num_particles(sp1), 2);

    BOOST_CHECK_EQUAL((*space).list_particles_within_radius(Real3(0.996, 0.5, 0.5), radius).size(), 2);
    (*space).remove_particle(pid1);
    {
        const std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
            retval = (*space).list_particles_within_radius(Real3(0.996, 0.5, 0.5), radius);
        BOOST_CHECK_EQUAL(retval.size(), 1);
        BOOST_CHECK_EQUAL(retval[0].first.first, pid2);
        BOOST_CHECK_CLOSE(retval[0].second, 0.002 - radius, 1e-6);
    }
}

BOOST_AUTO_TEST_CASE(ParticleSpace_test_list_particles_within_radius_boundary_xyz)
{
    std::unique_ptr<ParticleSpace> space(new particle_space_type(edge_lengths));
    SerialIDGenerator<ParticleID> pidgen;

    const ParticleID pid1 = pidgen();
    const ParticleID pid2 = pidgen();
    const Species sp1 = Species("A");

    BOOST_CHECK((*space).update_particle(pid1, Particle(sp1, Real3(0.0025, 0.0025, 0.0025), radius, 0)));
    BOOST_CHECK_EQUAL((*space).num_particles(sp1), 1);

    {
        const std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
            retval = (*space).list_particles_within_radius(Real3(0.9975, 0.9975, 0.9975), radius);
        BOOST_CHECK_EQUAL(retval.size(), 1);
        BOOST_CHECK_EQUAL(retval[0].first.first, pid1);
        BOOST_CHECK_CLOSE(retval[0].second, std::sqrt(3.0) * 0.005 - radius, 1e-6);
    }
    BOOST_CHECK_EQUAL((*space).list_particles_within_radius(Real3(0.9975, 0.9975, 0.9975), radius, pid1).size(), 0);
    BOOST_CHECK_EQUAL((*space).list_particles_within_radius(Real3(0.9965, 0.9965, 0.9965), radius).size(), 0);

    BOOST_CHECK((*space).update_particle(pid2, Particle(sp1, Real3(0.9965, 0.9965, 0.9965), radius, 0)));
    BOOST_CHECK_EQUAL((*space).num_particles(sp1), 2);

    BOOST_CHECK_EQUAL((*space).list_particles_within_radius(Real3(0.9975, 0.9975, 0.9975), radius).size(), 2);
    (*space).remove_particle(pid1);
    {
        const std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
            retval = (*space).list_particles_within_radius(Real3(0.9975, 0.9975, 0.9975), radius);
        BOOST_CHECK_EQUAL(retval.size(), 1);
        BOOST_CHECK_EQUAL(retval[0].first.first, pid2);
        BOOST_CHECK_CLOSE(retval[0].second, std::sqrt(3.0) * 0.001 - radius, 1e-6);
    }
}




BOOST_AUTO_TEST_CASE(ParticleSpaceRTreeImpl_test_constructor)
{
    std::unique_ptr<ParticleSpaceRTreeImpl> space(new ParticleSpaceRTreeImpl(edge_lengths));

    BOOST_CHECK_EQUAL(space->edge_lengths()[0], edge_lengths[0]);
    BOOST_CHECK_EQUAL(space->edge_lengths()[1], edge_lengths[1]);
    BOOST_CHECK_EQUAL(space->edge_lengths()[2], edge_lengths[2]);
}

BOOST_AUTO_TEST_SUITE_END()

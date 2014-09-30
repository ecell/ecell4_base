#define BOOST_TEST_MODULE "ParticleSpace_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/test_case_template.hpp>

#include <ecell4/core/types.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/ParticleSpace.hpp>
#include <ecell4/core/ParticleSpaceCellListImpl.hpp>

using namespace ecell4;

template<typename Timpl_>
void ParticleSpace_test_edge_lengths_template()
{
    const Real L(1e-6);
    Timpl_ target(Position3(L, L, L));

    const Position3 edge_lengths(target.edge_lengths());
    for (Position3::size_type dim(0); dim < 3; ++dim)
    {
        BOOST_CHECK_CLOSE(edge_lengths[dim], L, 1e-6);
    }
}

BOOST_AUTO_TEST_CASE(ParticleSpace_test_edge_lengths)
{
    ParticleSpace_test_edge_lengths_template<ParticleSpaceVectorImpl>();
    ParticleSpace_test_edge_lengths_template<ParticleSpaceCellListImpl>();
}

template<typename Timpl_>
void ParticleSpace_test_manipulate_particle_template()
{
    const Real L(1e-6);
    const Real L_2(L * 0.5);
    Timpl_ target(Position3(L, L, L));
    SerialIDGenerator<ParticleID> pidgen;

    const Species sp1("A"), sp2("B");
    const Real radius(2.5e-9), D(1e-12);
    const ParticleID pid(pidgen());
    BOOST_CHECK(
        target.update_particle(pid, Particle(sp1, Position3(0, 0, 0), radius, D)));
    BOOST_CHECK(target.has_particle(pid));
    BOOST_CHECK_EQUAL(target.num_particles(), 1);
    BOOST_CHECK_EQUAL(target.num_particles(sp1), 1);
    BOOST_CHECK_EQUAL(target.num_particles_exact(sp1), 1);
    BOOST_CHECK_EQUAL(target.num_particles(sp2), 0);
    BOOST_CHECK_EQUAL(target.num_particles_exact(sp2), 0);

    const std::pair<ParticleID, Particle>
        pid_particle_pair1(target.get_particle(pid));
    BOOST_CHECK_EQUAL(pid_particle_pair1.first, pid);
    BOOST_CHECK(pid_particle_pair1.second.species() == sp1);
    BOOST_CHECK_CLOSE(pid_particle_pair1.second.radius(), radius, 1e-6);
    BOOST_CHECK_CLOSE(pid_particle_pair1.second.D(), D, 1e-6);

    for (Position3::size_type dim(0); dim < 3; ++dim)
    {
        BOOST_CHECK_CLOSE(pid_particle_pair1.second.position()[dim], 0, 1e-6);
    }

    BOOST_CHECK(
        !target.update_particle(pid, Particle(sp2, Position3(L_2, L_2, L_2), radius, D)));
    BOOST_CHECK(target.has_particle(pid));
    BOOST_CHECK_EQUAL(target.num_particles(), 1);
    BOOST_CHECK_EQUAL(target.num_particles(sp1), 0);
    BOOST_CHECK_EQUAL(target.num_particles_exact(sp1), 0);
    BOOST_CHECK_EQUAL(target.num_particles(sp2), 1);
    BOOST_CHECK_EQUAL(target.num_particles_exact(sp2), 1);

    const std::pair<ParticleID, Particle>
        pid_particle_pair2(target.get_particle(pid));
    BOOST_CHECK_EQUAL(pid_particle_pair2.first, pid);
    BOOST_CHECK(pid_particle_pair2.second.species() == sp2);
    for (Position3::size_type dim(0); dim < 3; ++dim)
    {
        BOOST_CHECK_CLOSE(pid_particle_pair2.second.position()[dim], L_2, 1e-6);
    }

    target.remove_particle(pid);
    BOOST_CHECK(!target.has_particle(pid));
    BOOST_CHECK_EQUAL(target.num_particles(), 0);
    BOOST_CHECK_EQUAL(target.num_particles(sp2), 0);
    BOOST_CHECK_EQUAL(target.num_particles_exact(sp2), 0);
}

BOOST_AUTO_TEST_CASE(ParticleSpace_test_manupulate_particle)
{
    ParticleSpace_test_manipulate_particle_template<ParticleSpaceVectorImpl>();
    ParticleSpace_test_manipulate_particle_template<ParticleSpaceCellListImpl>();
}

template<typename Timpl_>
void ParticleSpace_test_manipulate_multiple_particles_template()
{
    const Real L(1e-6);
    Timpl_ target(Position3(L, L, L));
    SerialIDGenerator<ParticleID> pidgen;

    const Species sp1("A"), sp2("B");
    const Real radius(2.5e-9), D(1e-12);
    const ParticleID pid1(pidgen()), pid2(pidgen()), pid3(pidgen()),
        pid4(pidgen()), pid5(pidgen());
    // const Real L_6(L / 6);
    const Position3
        pos1(0, 0, 0),
        pos2(0, 0, 0),
        pos3(0, 0, 0),
        pos4(0, 0, 0),
        pos5(0, 0, 0);

    BOOST_CHECK(target.update_particle(pid1, Particle(sp1, pos1, radius, D)));
    BOOST_CHECK(target.update_particle(pid2, Particle(sp1, pos2, radius, D)));
    BOOST_CHECK(target.update_particle(pid3, Particle(sp2, pos3, radius, D)));
    BOOST_CHECK(target.update_particle(pid4, Particle(sp1, pos4, radius, D)));
    BOOST_CHECK(target.update_particle(pid5, Particle(sp2, pos5, radius, D)));

    BOOST_CHECK(target.has_particle(pid2));
    BOOST_CHECK(target.has_particle(pid3));
    BOOST_CHECK_EQUAL(target.num_particles(), 5);
    BOOST_CHECK_EQUAL(target.num_particles(sp1), 3);
    BOOST_CHECK_EQUAL(target.num_particles(sp2), 2);
    BOOST_CHECK_EQUAL(target.num_particles_exact(sp1), 3);
    BOOST_CHECK_EQUAL(target.num_particles_exact(sp2), 2);
    BOOST_CHECK_EQUAL(target.list_particles().size(), 5);
    BOOST_CHECK_EQUAL(target.list_particles(sp1).size(), 3);
    BOOST_CHECK_EQUAL(target.list_particles(sp2).size(), 2);
    BOOST_CHECK_EQUAL(target.list_particles_exact(sp1).size(), 3);
    BOOST_CHECK_EQUAL(target.list_particles_exact(sp2).size(), 2);

    target.remove_particle(pid2);
    target.remove_particle(pid3);

    BOOST_CHECK(!target.has_particle(pid2));
    BOOST_CHECK(!target.has_particle(pid3));
    BOOST_CHECK_EQUAL(target.num_particles(), 3);
    BOOST_CHECK_EQUAL(target.num_particles(sp1), 2);
    BOOST_CHECK_EQUAL(target.num_particles(sp2), 1);
}

BOOST_AUTO_TEST_CASE(ParticleSpace_test_multiple_particles)
{
    ParticleSpace_test_manipulate_multiple_particles_template<ParticleSpaceVectorImpl>();
    ParticleSpace_test_manipulate_multiple_particles_template<ParticleSpaceCellListImpl>();
}

template<typename Timpl_>
void ParticleSpace_test_list_particles_within_radius_template()
{
    const Real L(1e-6);
    Timpl_ target(Position3(L, L, L));
    SerialIDGenerator<ParticleID> pidgen;

    const Species sp1("A"), sp2("B");
    const Real radius(2.5e-9), D(1e-12);
    const ParticleID pid1(pidgen()), pid2(pidgen()), pid3(pidgen()),
        pid4(pidgen()), pid5(pidgen()), pid6(pidgen());
    const Real L_8(L * 0.125), L_4(L * 0.25);
    const Position3
        pos1(L_8, L_8, L_8),
        // pos2(L * -0.04, L * -0.04, L * -0.04),
        pos2(L * 0.96, L * 0.96, L * 0.96),
        pos3(L_8 + L_4, L_8, L_8),
        pos4(L_8 + L_4, L_8 + L_4, L_8),
        pos5(L_4, L_4, L_4),
        pos6(L_4 + L_8, L_4 + L_8, L_4 + L_8);

    BOOST_CHECK(target.update_particle(pid1, Particle(sp1, pos1, radius, D)));
    BOOST_CHECK(target.update_particle(pid2, Particle(sp1, pos2, radius, D)));
    BOOST_CHECK(target.update_particle(pid3, Particle(sp2, pos3, radius, D)));
    BOOST_CHECK(target.update_particle(pid4, Particle(sp1, pos4, radius, D)));
    BOOST_CHECK(target.update_particle(pid5, Particle(sp2, pos5, radius, D)));
    BOOST_CHECK(target.update_particle(pid6, Particle(sp2, pos6, radius, D)));

    BOOST_CHECK_EQUAL(target.num_particles(), 6);

    const Real R1(L * 0.3);
    const std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        retval1(target.list_particles_within_radius(pos1, R1));
    BOOST_CHECK_EQUAL(retval1.size(), 4);
    BOOST_CHECK_CLOSE(retval1[0].second, 0, 1e-6);
    BOOST_CHECK_EQUAL(retval1[0].first.first, pid1);
    BOOST_CHECK_EQUAL(retval1[1].first.first, pid5);
    BOOST_CHECK_EQUAL(retval1[2].first.first, pid3);
    BOOST_CHECK_EQUAL(retval1[3].first.first, pid2);
    BOOST_CHECK_CLOSE(retval1[3].second, 0.2857883832488648 * L, 1e-6);

    BOOST_CHECK_EQUAL(target.list_particles_within_radius(pos1, R1, pid1).size(), 3);
    BOOST_CHECK_EQUAL(target.list_particles_within_radius(pos1, R1, pid1, pid3).size(), 2);

    const Real R2(L * 0.3);
    const std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        retval2(target.list_particles_within_radius(pos5, R2, pid5));
    BOOST_CHECK_EQUAL(retval2.size(), 4);
    BOOST_CHECK(
        retval2[0].first.first == pid1 || retval2[0].first.first == pid6
        || retval2[0].first.first == pid4 || retval2[0].first.first == pid3);
    BOOST_CHECK(
        retval2[1].first.first == pid1 || retval2[1].first.first == pid6
        || retval2[1].first.first == pid4 || retval2[1].first.first == pid3);
    BOOST_CHECK(
        retval2[2].first.first == pid1 || retval2[2].first.first == pid6
        || retval2[2].first.first == pid4 || retval2[2].first.first == pid3);
    BOOST_CHECK(
        retval2[3].first.first == pid1 || retval2[3].first.first == pid6
        || retval2[3].first.first == pid4 || retval2[3].first.first == pid3);
    BOOST_CHECK_CLOSE(retval2[0].second, retval2[1].second, 1e-6);
    BOOST_CHECK_CLOSE(retval2[0].second, retval2[2].second, 1e-6);
    BOOST_CHECK_CLOSE(retval2[0].second, retval2[3].second, 1e-6);
}

BOOST_AUTO_TEST_CASE(ParticleSpace_test_list_particles_within_radius)
{
    ParticleSpace_test_list_particles_within_radius_template<ParticleSpaceVectorImpl>();
    ParticleSpace_test_list_particles_within_radius_template<ParticleSpaceCellListImpl>();
}

#define BOOST_TEST_MODULE "OffLatticeSpace_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/floating_point_comparison.hpp>

#include <ecell4/core/OffLatticeSpace.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>

using namespace ecell4;

struct Fixture
{
    const Real voxel_radius;
    const Species species;
    const Voxel voxel;
    OffLatticeSpace space;
    SerialIDGenerator<ParticleID> sidgen;

    Fixture() :
        voxel_radius(2.5e-9),
        species(/* serial = */ "SpeciesA",
                /* radius = */ "2.5e-9",
                /* D = */      "1e-12"),
        voxel(/* species = */    species,
              /* coordinate = */ 3,
              /* radius = */     2.5e-9,
              /* D = */          1e-12),
        space(voxel_radius)
    {
        OffLatticeSpace::position_container positions;
        const Real unit(voxel_radius / sqrt(3.0));
        for (int i(0); i < 10; ++i)
            positions.push_back(
                    Real3(unit * i, unit * i, unit * i));
        OffLatticeSpace::coordinate_pair_list_type adjoining_pairs;
        for (int i(1); i < 10; ++i )
            adjoining_pairs.push_back(
                    std::make_pair(i-1, i));
        space = OffLatticeSpace(voxel_radius, positions, adjoining_pairs);
    }
};

BOOST_FIXTURE_TEST_SUITE(suite, Fixture)

BOOST_AUTO_TEST_CASE(OffLatticeSpace_test_constructor) {}

BOOST_AUTO_TEST_CASE(OffLatticeSpace_test_molecules)
{
    const ParticleID pid(sidgen());
    space.update_voxel(pid, voxel);

    BOOST_CHECK_EQUAL(space.num_molecules(species), 1);
}

BOOST_AUTO_TEST_CASE(OffLatticeSpace_test_voxelspacebase)
{
    const ParticleID pid(sidgen());
    space.update_voxel(pid, voxel);

    BOOST_CHECK_EQUAL(space.list_species().size(), 1);
    BOOST_CHECK_EQUAL(space.num_voxels_exact(species), 1);
    BOOST_CHECK_EQUAL(space.num_voxels(species), 1);
    BOOST_CHECK_EQUAL(space.num_voxels(), 1);

    BOOST_CHECK(space.has_voxel(pid));
    BOOST_CHECK(!space.has_voxel(sidgen()));

    BOOST_CHECK_EQUAL(space.list_voxels().size(), 1);
    BOOST_CHECK_EQUAL(space.list_voxels(species).size(), 1);
    BOOST_CHECK_EQUAL(space.list_voxels_exact(species).size(), 1);

    BOOST_CHECK_EQUAL(space.list_voxels().at(0).first, pid);
    BOOST_CHECK_EQUAL(space.list_voxels(species).at(0).first, pid);
    BOOST_CHECK_EQUAL(space.list_voxels_exact(species).at(0).first, pid);

    BOOST_CHECK_EQUAL(space.get_voxel(pid).first, pid);

    BOOST_CHECK_NO_THROW(space.find_voxel_pool(species));
    BOOST_CHECK(space.has_molecule_pool(species));
    BOOST_CHECK_NO_THROW(space.find_molecule_pool(species));
}

BOOST_AUTO_TEST_CASE(OffLatticeSpace_test_voxel)
{
    const ParticleID pid(sidgen());

    BOOST_CHECK(space.update_voxel(pid, voxel));
    BOOST_CHECK(space.remove_voxel(pid));
    BOOST_CHECK(!space.remove_voxel(3));

    BOOST_CHECK(space.update_voxel(pid, voxel));
    BOOST_CHECK(space.remove_voxel(3));
    BOOST_CHECK(!space.remove_voxel(pid));
}

BOOST_AUTO_TEST_CASE(OffLatticeSpace_test_move)
{
    const ParticleID pid(sidgen());
    BOOST_CHECK(space.update_voxel(pid, voxel));

    BOOST_CHECK(space.can_move(3, 4));
    BOOST_CHECK(space.move(3, 4, 0));
    BOOST_CHECK_EQUAL(space.get_voxel_at(3).first, ParticleID());
    BOOST_CHECK_EQUAL(space.get_voxel_at(4).first, pid);

    BOOST_CHECK(!space.can_move(3, 4));
}

BOOST_AUTO_TEST_CASE(OffLatticeSpace_test_at)
{
    BOOST_CHECK_EQUAL(space.size(), 10);

    const ParticleID pid(sidgen());
    BOOST_CHECK(space.update_voxel(pid, voxel));

    BOOST_CHECK_NO_THROW(space.get_voxel_at(3));
    BOOST_CHECK_EQUAL(space.get_voxel_at(3).first, pid);

    BOOST_CHECK_NO_THROW(space.particle_at(3));
    BOOST_CHECK_EQUAL(space.particle_at(3).species(), species);
    BOOST_CHECK_EQUAL(space.particle_at(3).position(), space.coordinate2position(3));
    BOOST_CHECK_EQUAL(space.particle_at(3).radius(), 2.5e-9);
    BOOST_CHECK_EQUAL(space.particle_at(3).D(), 1e-12);
}

BOOST_AUTO_TEST_CASE(OffLatticeSpace_test_neighbor)
{
    BOOST_CHECK_EQUAL(space.num_neighbors(0), 1);
    BOOST_CHECK_EQUAL(space.num_neighbors(1), 2);
    BOOST_CHECK_EQUAL(space.num_neighbors(2), 2);
    BOOST_CHECK_EQUAL(space.num_neighbors(3), 2);
    BOOST_CHECK_EQUAL(space.num_neighbors(4), 2);
    BOOST_CHECK_EQUAL(space.num_neighbors(5), 2);
    BOOST_CHECK_EQUAL(space.num_neighbors(6), 2);
    BOOST_CHECK_EQUAL(space.num_neighbors(7), 2);
    BOOST_CHECK_EQUAL(space.num_neighbors(8), 2);
    BOOST_CHECK_EQUAL(space.num_neighbors(9), 1);

    BOOST_CHECK_EQUAL(space.get_neighbor(0, 0), 1);
    BOOST_CHECK_EQUAL(space.get_neighbor(1, 0), 0);
    BOOST_CHECK_EQUAL(space.get_neighbor(1, 1), 2);
    BOOST_CHECK_EQUAL(space.get_neighbor(2, 0), 1);
    BOOST_CHECK_EQUAL(space.get_neighbor(2, 1), 3);
    BOOST_CHECK_EQUAL(space.get_neighbor(3, 0), 2);
    BOOST_CHECK_EQUAL(space.get_neighbor(3, 1), 4);
    BOOST_CHECK_EQUAL(space.get_neighbor(4, 0), 3);
    BOOST_CHECK_EQUAL(space.get_neighbor(4, 1), 5);
    BOOST_CHECK_EQUAL(space.get_neighbor(5, 0), 4);
    BOOST_CHECK_EQUAL(space.get_neighbor(5, 1), 6);
    BOOST_CHECK_EQUAL(space.get_neighbor(6, 0), 5);
    BOOST_CHECK_EQUAL(space.get_neighbor(6, 1), 7);
    BOOST_CHECK_EQUAL(space.get_neighbor(7, 0), 6);
    BOOST_CHECK_EQUAL(space.get_neighbor(7, 1), 8);
    BOOST_CHECK_EQUAL(space.get_neighbor(8, 0), 7);
    BOOST_CHECK_EQUAL(space.get_neighbor(8, 1), 9);
    BOOST_CHECK_EQUAL(space.get_neighbor(9, 0), 8);
}

BOOST_AUTO_TEST_SUITE_END()

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
    OffLatticeSpace space;
    SerialIDGenerator<ParticleID> sidgen;

    Fixture() :
        voxel_radius(2.5e-9), space(voxel_radius)
    {
        OffLatticeSpace::position_container positions;
        const Real unit(voxel_radius / sqrt(3.0));
        for (int i(0); i < 10; ++i)
            positions.push_back(Real3(unit * i, unit * i, unit * i));
        OffLatticeSpace::coordinate_pair_list_type adjoining_pairs;
        for (int i(1); i < 10; ++i )
            adjoining_pairs.push_back(std::make_pair(i-1, i));
        space = OffLatticeSpace(voxel_radius, positions, adjoining_pairs);
    }
};

BOOST_FIXTURE_TEST_SUITE(suite, Fixture)

BOOST_AUTO_TEST_CASE(OffLatticeSpace_test_constructor) {}

BOOST_AUTO_TEST_CASE(OffLatticeSpace_test_something)
{
    BOOST_CHECK_EQUAL(space.size(), 10);

    const Species species("SpeciesA", "2.5e-9", "1e-12");
    const Voxel voxel(species, 3, /* r= */2.5e-9, /* D= */1e-12);
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

    BOOST_CHECK_NO_THROW(space.find_molecule_pool(species));
}

BOOST_AUTO_TEST_SUITE_END()

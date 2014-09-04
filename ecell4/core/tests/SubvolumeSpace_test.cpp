#define BOOST_TEST_MODULE "SubvolumeSpace_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/test_case_template.hpp>

#include <ecell4/core/types.hpp>
#include <ecell4/core/SubvolumeSpace.hpp>

using namespace ecell4;

template<typename Timpl_>
void SubvolumeSpace_test_volume_template()
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Global matrix_sizes(2, 3, 4);
    Timpl_ target(edge_lengths, matrix_sizes);

    BOOST_CHECK_CLOSE(target.volume(), 1e-18, 1e-6);
    BOOST_CHECK_EQUAL(target.num_subvolumes(), 24);
    BOOST_CHECK_CLOSE(target.subvolume(), 4.166666666666667e-20, 1e-6);

    BOOST_CHECK_EQUAL(target.global2coord(Global()), 0);
    BOOST_CHECK_EQUAL(target.global2coord(Global(1, 0, 0)), 1);
    BOOST_CHECK_EQUAL(target.global2coord(Global(1, 2, 3)), 23);
    BOOST_CHECK_EQUAL(target.coord2global(0), Global(0, 0, 0));
    BOOST_CHECK_EQUAL(target.coord2global(1), Global(1, 0, 0));
    BOOST_CHECK_EQUAL(target.coord2global(23), Global(1, 2, 3));
}

BOOST_AUTO_TEST_CASE(SubvolumeSpace_test_volume)
{
    SubvolumeSpace_test_volume_template<SubvolumeSpaceVectorImpl>();
}

template<typename Timpl_>
void SubvolumeSpace_test_num_molecules_template()
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Global matrix_sizes(2, 3, 4);
    Timpl_ target(edge_lengths, matrix_sizes);

    const Species sp1("A"), sp2("B");
    BOOST_CHECK_EQUAL(target.num_molecules_exact(sp1, 0), 0);
    BOOST_CHECK_EQUAL(target.num_molecules_exact(sp1, 23), 0);
    target.add_molecules(sp1, 60, 0);
    target.remove_molecules(sp1, 30, 0);
    target.add_molecules(sp1, 60, 23);
    target.add_molecules(sp2, 60, 23);
    BOOST_CHECK_EQUAL(target.num_molecules_exact(sp1, 0), 30);
    BOOST_CHECK_EQUAL(target.num_molecules_exact(sp2, 0), 0);
    BOOST_CHECK_EQUAL(target.num_molecules_exact(sp1, 23), 60);
    BOOST_CHECK_EQUAL(target.num_molecules_exact(sp2, 23), 60);

    BOOST_CHECK_EQUAL(target.num_molecules_exact(sp1), 90);
    BOOST_CHECK_EQUAL(target.num_molecules_exact(sp2), 60);
    BOOST_CHECK_EQUAL(target.num_molecules(Species("_")), 150);
}

BOOST_AUTO_TEST_CASE(SubvolumeSpace_test_num_molecules)
{
    SubvolumeSpace_test_num_molecules_template<SubvolumeSpaceVectorImpl>();
}

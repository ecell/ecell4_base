#define BOOST_TEST_MODULE "CompartmentSpace_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>

#include "../types.hpp"
#include "../CompartmentSpace.hpp"

using namespace ecell4;

template<typename Timpl_>
void CompartmentSpace_test_volume_template()
{
    const Real volume(1e-18);
    Timpl_ target(volume);
    const Real new_volume(2 * target.volume());
    target.set_volume(new_volume);
    BOOST_CHECK_EQUAL(target.volume(), new_volume);
}

BOOST_AUTO_TEST_CASE(CompartmentSpace_test_volume)
{
    CompartmentSpace_test_volume_template<CompartmentSpaceVectorImpl>();
}

template<typename Timpl_>
void CompartmentSpace_test_species_template()
{
    const Real volume(1e-18);
    Timpl_ target(volume);

    Species sp1("A"), sp2("B"), sp3("C");
    target.add_species(sp1);
    target.add_species(sp2);
    BOOST_CHECK(target.has_species(sp1));
    BOOST_CHECK(target.has_species(sp2));
    BOOST_CHECK(!target.has_species(sp3));
    target.add_species(sp3);
    BOOST_CHECK(target.has_species(sp3));
    BOOST_CHECK(target.num_species() == 3);

    target.remove_species(sp2);
    BOOST_CHECK(!target.has_species(sp2));
    BOOST_CHECK_THROW(target.remove_species(sp2), NotFound);
    BOOST_CHECK(target.has_species(sp1));
    BOOST_CHECK(target.has_species(sp3));

    BOOST_CHECK(target.num_molecules(sp1) == 0);
    target.add_molecules(sp1, 30);
    BOOST_CHECK(target.num_molecules(sp1) == 30);
    target.remove_molecules(sp1, 10);
    BOOST_CHECK(target.num_molecules(sp1) == 20);
}

BOOST_AUTO_TEST_CASE(CompartmentSpace_test_species)
{
    CompartmentSpace_test_species_template<CompartmentSpaceVectorImpl>();
}

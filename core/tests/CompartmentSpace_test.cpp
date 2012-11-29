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
    Real const volume(1e-18);
    Timpl_ target(volume);
    target.set_volume(2 * target.volume());
}

BOOST_AUTO_TEST_CASE(CompartmentSpace_test_volume)
{
    CompartmentSpace_test_volume_template<CompartmentSpaceVectorImpl>();
}


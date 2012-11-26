#define BOOST_TEST_MODULE "CompartmentSpace_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "../Position3.hpp"
#include "../linear_algebra.hpp"
#include "../types.hpp"
#include "../CompartmentSpace.hpp"

using namespace ecell4;


//____________________________________________________________________________//

typedef boost::mpl::list<CompartmentSpaceVectorImpl> test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(CompartmentSpace_test_volume1, Timpl_, test_types)
{
    Real const volume(1e-18);
    Timpl_ target(volume);
    target.set_volume(2 * target.volume());
}

template<typename Timpl_>
void CompartmentSpace_test_volume2()
{
    Real const volume(1e-18);
    Timpl_ target(volume);
    target.set_volume(2 * target.volume());
}

BOOST_AUTO_TEST_CASE(CompartmentSpace_test_volume)
{
    CompartmentSpace_test_volume2<CompartmentSpaceVectorImpl>();
}

//____________________________________________________________________________//

BOOST_AUTO_TEST_CASE( Position3Multiply )
{
  Position3 pos1(1,2,3);
  BOOST_CHECK_EQUAL( pos1 * 2, Position3(2,4,6));
}

BOOST_AUTO_TEST_CASE( Position3Add )
{
  Position3 pos2(1,2,3);
  Position3 pos3(2,4,6);
  BOOST_CHECK_EQUAL( pos2 + pos3, Position3(3,6,9));
}

BOOST_AUTO_TEST_CASE( Position3Sub )
{
  Position3 pos4(2,4,6);
  Position3 pos5(1,2,3);
  BOOST_CHECK_EQUAL( pos4 - pos5, Position3(1,2,3));
}

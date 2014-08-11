#define BOOST_TEST_MODULE "EventScheduler_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/floating_point_comparison.hpp>

#include <ecell4/core/EventScheduler.hpp>

using namespace ecell4;

BOOST_AUTO_TEST_CASE(EventScheduler_test_constructor)
{
    EventScheduler scheduler;
}

#define BOOST_TEST_MODULE "extras_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <ecell4/core/extras.hpp>

using namespace ecell4;


BOOST_AUTO_TEST_CASE(extras_test_)
{
    const extras::VersionInformation vinfo = extras::parse_version_information("ecell4-test-1.2.3");

    BOOST_CHECK_EQUAL(vinfo.header, "ecell4-test-");
    BOOST_CHECK_EQUAL(vinfo.majorno, 1);
    BOOST_CHECK_EQUAL(vinfo.minorno, 2);
    BOOST_CHECK_EQUAL(vinfo.patchno, 3);
}

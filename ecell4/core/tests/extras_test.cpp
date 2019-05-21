#define BOOST_TEST_MODULE "extras_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <ecell4/core/extras.hpp>
#include <ecell4/core/NetworkModel.hpp>

using namespace ecell4;


BOOST_AUTO_TEST_CASE(extras_test_)
{
    const extras::VersionInformation vinfo = extras::parse_version_information("ecell4-test-1.2.3");

    BOOST_CHECK_EQUAL(vinfo.header, "ecell4-test-");
    BOOST_CHECK_EQUAL(vinfo.majorno, 1);
    BOOST_CHECK_EQUAL(vinfo.minorno, 2);
    BOOST_CHECK_EQUAL(vinfo.patchno, 3);
}

BOOST_AUTO_TEST_CASE(DimensionAttributeTest)
{
    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_species_attribute(Species("A"));
    model->add_species_attribute(Species("B", 1.0, 0.0, "A"));
    model->add_species_attribute(Species("C", 1.0, 0.0, "", 2));
    model->add_species_attribute(Species("D", 1.0, 0.0, "C"));

    BOOST_CHECK_EQUAL(extras::get_dimension_from_model(Species("A"), model), Shape::THREE);
    BOOST_CHECK_EQUAL(extras::get_dimension_from_model(Species("B"), model), Shape::THREE);
    BOOST_CHECK_EQUAL(extras::get_dimension_from_model(Species("C"), model), Shape::TWO);
    BOOST_CHECK_EQUAL(extras::get_dimension_from_model(Species("D"), model), Shape::TWO);
    BOOST_CHECK_THROW(extras::get_dimension_from_model(Species("E"), model), NotFound);
}

BOOST_AUTO_TEST_CASE(VersionInformationTest)
{
    {
        const extras::VersionInformation vinfo1 = extras::parse_version_information("ecell4-test-1.2.3");
        BOOST_CHECK_EQUAL(vinfo1.header, "ecell4-test-");
        BOOST_CHECK_EQUAL(vinfo1.majorno, 1);
        BOOST_CHECK_EQUAL(vinfo1.minorno, 2);
        BOOST_CHECK_EQUAL(vinfo1.patchno, 3);
        BOOST_CHECK_EQUAL(vinfo1.devno, -1);
    }
    {
        const extras::VersionInformation vinfo1 = extras::parse_version_information("ecell4-test-1.2");
        BOOST_CHECK_EQUAL(vinfo1.header, "ecell4-test-");
        BOOST_CHECK_EQUAL(vinfo1.majorno, 1);
        BOOST_CHECK_EQUAL(vinfo1.minorno, 2);
        BOOST_CHECK_EQUAL(vinfo1.patchno, 0);
        BOOST_CHECK_EQUAL(vinfo1.devno, -1);
    }
    {
        const extras::VersionInformation vinfo1 = extras::parse_version_information("ecell4-test-1.2.dev4");
        BOOST_CHECK_EQUAL(vinfo1.header, "ecell4-test-");
        BOOST_CHECK_EQUAL(vinfo1.majorno, 1);
        BOOST_CHECK_EQUAL(vinfo1.minorno, 2);
        BOOST_CHECK_EQUAL(vinfo1.patchno, 0);
        BOOST_CHECK_EQUAL(vinfo1.devno, 4);
    }
    {
        const extras::VersionInformation vinfo1 = extras::parse_version_information("ecell4-test-1.2.3c4.dev5");
        BOOST_CHECK_EQUAL(vinfo1.header, "ecell4-test-");
        BOOST_CHECK_EQUAL(vinfo1.majorno, 1);
        BOOST_CHECK_EQUAL(vinfo1.minorno, 2);
        BOOST_CHECK_EQUAL(vinfo1.patchno, 3);
        BOOST_CHECK_EQUAL(vinfo1.pre, extras::VersionInformation::RC);
        BOOST_CHECK_EQUAL(vinfo1.preno, 4);
        BOOST_CHECK_EQUAL(vinfo1.devno, 5);
    }
    {
        BOOST_CHECK(extras::check_version_information("ecell4-test-1.0", "ecell4-test-1.0"));
        BOOST_CHECK(extras::check_version_information("ecell4-test-1.1", "ecell4-test-1.0"));
        BOOST_CHECK(extras::check_version_information("ecell4-test-2.0.0", "ecell4-test-1.0"));
        BOOST_CHECK(extras::check_version_information("ecell4-test-2.0.0b1", "ecell4-test-1.0"));
        BOOST_CHECK(extras::check_version_information("ecell4-test-2.0.0b1", "ecell4-test-2.0.0b1"));
        BOOST_CHECK(extras::check_version_information("ecell4-test-1.0.1", "ecell4-test-2.0.0b1"));
        BOOST_CHECK(extras::check_version_information("ecell4-test-1.0", "ecell4-test-2.0.0b1"));
        BOOST_CHECK(!extras::check_version_information("ecell4-test-1.0.dev1", "ecell4-test-1.0"));
        BOOST_CHECK(extras::check_version_information("ecell4-test-1.1.dev1", "ecell4-test-1.0"));
        BOOST_CHECK(extras::check_version_information("ecell4-test-1.0.dev1", "ecell4-test-1.0.dev1"));
        BOOST_CHECK(extras::check_version_information("ecell4-test-1.0.dev2", "ecell4-test-1.0.dev1"));
        BOOST_CHECK(extras::check_version_information("ecell4-test-2.0.0b1.dev1", "ecell4-test-1.0.dev1"));
        BOOST_CHECK(!extras::check_version_information("ecell4-test-1.0a1", "ecell4-test-1.0"));
        BOOST_CHECK(extras::check_version_information("ecell4-test-1.1a2", "ecell4-test-1.0"));
        BOOST_CHECK(extras::check_version_information("ecell4-test-1.0rc1", "ecell4-test-1.0a5"));
        BOOST_CHECK(extras::check_version_information("ecell4-test-1.0.1a1", "ecell4-test-1.0"));
        BOOST_CHECK(extras::check_version_information("ecell4-test-1.0.1c1", "ecell4-test-1.0.1rc1"));
    }
}

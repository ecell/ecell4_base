#define BOOST_TEST_MODULE "Species_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <iostream>

#include <ecell4/core/Species.hpp>
#include <ecell4/core/Context.hpp>

using namespace ecell4;


BOOST_AUTO_TEST_CASE(Species_test_constructor)
{
    Species sp1("A");
    Species sp2("A(  )");
//     Species sp2("A( a   , b ^ 1,   c)");
    Species sp3("A(a,b^1,c=p).B(a,b=u^1)");
}

BOOST_AUTO_TEST_CASE(Species_test_name)
{
    Species sp1("test");
    BOOST_CHECK_EQUAL(sp1.serial(), "test");
    Species sp2("test()");
    BOOST_CHECK_EQUAL(sp2.serial(), "test()");
    // BOOST_CHECK_EQUAL(sp2.serial(), "test");
}

BOOST_AUTO_TEST_CASE(Species_test_attributes)
{
    Species sp("test");

    const std::string val1 = "foo";
    const Real val2 = 1.5;
    const Integer val3 = 3;
    const bool val4 = true;

    sp.set_attribute("attr1", val1);
    sp.set_attribute("attr2", val2);
    sp.set_attribute("attr3", val3);
    sp.set_attribute("attr4", val4);
    sp.set_attribute("attr5", "bar");

    BOOST_ASSERT(sp.has_attribute("attr1"));
    BOOST_ASSERT(sp.has_attribute("attr2"));
    BOOST_ASSERT(sp.has_attribute("attr3"));
    BOOST_ASSERT(sp.has_attribute("attr4"));
    BOOST_ASSERT(sp.has_attribute("attr5"));
    BOOST_ASSERT(!sp.has_attribute("atrt1"));

    const Species::attribute_type attr1 = sp.get_attribute("attr1");
    const Species::attribute_type attr2 = sp.get_attribute("attr2");
    const Species::attribute_type attr3 = sp.get_attribute("attr3");
    const Species::attribute_type attr4 = sp.get_attribute("attr4");
    const Species::attribute_type attr5 = sp.get_attribute("attr5");

    BOOST_ASSERT(attr1.type() == typeid(std::string));
    BOOST_ASSERT(attr2.type() == typeid(Real));
    BOOST_ASSERT(attr3.type() == typeid(Integer));
    BOOST_ASSERT(attr4.type() == typeid(bool));
    BOOST_ASSERT(attr5.type() == typeid(std::string));

    BOOST_CHECK_EQUAL(boost::get<std::string>(attr1), val1);
    BOOST_CHECK_EQUAL(boost::get<Real>(attr2), val2);
    BOOST_CHECK_EQUAL(boost::get<Integer>(attr3), val3);
    BOOST_CHECK_EQUAL(boost::get<bool>(attr4), val4);
    BOOST_CHECK_EQUAL(boost::get<std::string>(attr5), "bar");

    BOOST_CHECK_EQUAL(sp.get_attribute_as<std::string>("attr1"), val1);
    BOOST_CHECK_EQUAL(sp.get_attribute_as<Real>("attr2"), val2);
    BOOST_CHECK_EQUAL(sp.get_attribute_as<Integer>("attr3"), val3);
    BOOST_CHECK_EQUAL(sp.get_attribute_as<bool>("attr4"), val4);
    BOOST_CHECK_EQUAL(sp.get_attribute_as<std::string>("attr5"), "bar");

    sp.remove_attribute("attr1");
    sp.remove_attribute("attr2");
    sp.remove_attribute("attr3");
    sp.remove_attribute("attr4");
    sp.remove_attribute("attr5");

    BOOST_ASSERT(!sp.has_attribute("attr1"));
    BOOST_ASSERT(!sp.has_attribute("attr2"));
    BOOST_ASSERT(!sp.has_attribute("attr3"));
    BOOST_ASSERT(!sp.has_attribute("attr4"));
    BOOST_ASSERT(!sp.has_attribute("attr5"));

    // Species species("test");
    // species.set_attribute("attr1", "value1");
    // species.set_attribute("attr2", "value2");
    // BOOST_CHECK_EQUAL(species.get_attribute_as<std::string>("attr1"), "value1");
    // BOOST_CHECK_EQUAL(species.get_attribute_as<std::string>("attr2"), "value2");
    // species.remove_attribute("attr1");
}

BOOST_AUTO_TEST_CASE(Species_test_match1)
{
    Species sp1, sp2;

    sp1.add_unit(UnitSpecies("C"));
    sp1.add_unit(UnitSpecies("A"));
    sp1.add_unit(UnitSpecies("B"));
    BOOST_CHECK_EQUAL(sp1.serial(), "C.A.B");

    sp2.add_unit(UnitSpecies("A"));
    sp2.add_unit(UnitSpecies("C"));
    BOOST_CHECK_EQUAL(sp2.serial(), "A.C");

    // BOOST_CHECK(sp2.match(sp1));
    // BOOST_CHECK(!sp1.match(sp2));
}

BOOST_AUTO_TEST_CASE(Species_test_match2)
{
    Species sp1, sp2;

    sp1.add_unit(UnitSpecies("B"));
    sp1.add_unit(UnitSpecies("A"));
    sp1.add_unit(UnitSpecies("A"));
    BOOST_CHECK_EQUAL(sp1.serial(), "B.A.A");

    sp2.add_unit(UnitSpecies("B"));
    sp2.add_unit(UnitSpecies("A"));
    BOOST_CHECK_EQUAL(sp2.serial(), "B.A");

    // BOOST_CHECK(sp2.match(sp1));
    // BOOST_CHECK(!sp1.match(sp2));
}

BOOST_AUTO_TEST_CASE(Species_test_serialization)
{
    Species sp1("X(a^1).Y(a^3,b).X(a^2).Y(a^1,b^2).X(a^3)");

    BOOST_CHECK_EQUAL(
        format_species(Species("X(a^1).Y(a^3,b).X(a^2).Y(a^1,b^2).X(a^3)")).serial(),
        "X(a^1).Y(a^1,b).X(a^2).Y(a^2,b^3).X(a^3)");
    BOOST_CHECK_EQUAL(
        format_species(Species("X(a^1).Y(a^3,b^4).X(a^3).Z(a^4,b^5).Y(a^1,b^2).Z(a^2,b^5)")).serial(),
        "X(a^1).Y(a^1,b^2).Z(a^2,b^3).Z(a^4,b^3).Y(a^5,b^4).X(a^5)");
    BOOST_CHECK_EQUAL(
        format_species(Species("X(a^3,b^1).X(a^2,b).X(a,b^3).X(a^1,b^4).X(a^4,b^2)")).serial(),
        "X(a,b^1).X(a^1,b^2).X(a^2,b^3).X(a^3,b^4).X(a^4,b)");
}

BOOST_AUTO_TEST_CASE(Species_test_match3)
{
    BOOST_CHECK_EQUAL(count_spmatches(Species("A"), Species("A.A.A")), 3);

    BOOST_CHECK_EQUAL(count_spmatches(Species("_1._2"), Species("A.B.C")), 6);

    MatchObject::context_type::variable_container_type globals;
    globals["_1"] = "A";
    BOOST_CHECK_EQUAL(
        count_spmatches(Species("_1._2"), Species("A.B.C"), globals), 2);
}

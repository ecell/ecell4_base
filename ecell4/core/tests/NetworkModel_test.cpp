#define BOOST_TEST_MODULE "NetworkModel_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/NetworkModel.hpp>

using namespace ecell4;


BOOST_AUTO_TEST_CASE(NetworkModel_test_constructor)
{
    NetworkModel model;
}

BOOST_AUTO_TEST_CASE(NetworkModel_test_species)
{
    NetworkModel model;

    {
        Species sp1("A");
        sp1.set_attribute("key1", "value1");
        Species sp2("B");
        sp2.set_attribute("key1", "value2");
        Species sp3("_");
        sp3.set_attribute("key1", "value0");

        model.add_species_attribute(sp1);
        model.add_species_attribute(sp2);
        model.add_species_attribute(sp3);

        BOOST_CHECK(model.has_species_attribute(Species("_")));
        BOOST_CHECK(model.has_species_attribute(Species("A")));
        BOOST_CHECK(model.has_species_attribute(Species("B")));
        BOOST_CHECK(!model.has_species_attribute(Species("C")));

        Species sp4("A");
        sp4.set_attribute("key2", "value3");
        BOOST_CHECK(!model.update_species_attribute(sp4));
    }

    {
        const Species sp1 = model.apply_species_attributes(Species("A"));
        BOOST_CHECK(sp1.has_attribute("key1"));
        BOOST_CHECK_EQUAL(sp1.get_attribute_as<std::string>("key1"), "value1");
        BOOST_CHECK(sp1.has_attribute("key2"));
        BOOST_CHECK_EQUAL(sp1.get_attribute_as<std::string>("key2"), "value3");

        const Species sp2 = model.apply_species_attributes(Species("B"));
        BOOST_CHECK(sp2.has_attribute("key1"));
        BOOST_CHECK_EQUAL(sp2.get_attribute_as<std::string>("key1"), "value2");
        BOOST_CHECK(!sp2.has_attribute("key2"));

        const Species sp3 = model.apply_species_attributes(Species("C"));
        BOOST_CHECK(sp3.has_attribute("key1"));
        BOOST_CHECK_EQUAL(sp3.get_attribute_as<std::string>("key1"), "value0");
        BOOST_CHECK(!sp3.has_attribute("key2"));
    }

    model.remove_species_attribute(Species("A"));
    BOOST_CHECK(!model.has_species_attribute(Species("A")));
    BOOST_CHECK_THROW(model.remove_species_attribute(Species("A")), NotFound);
}

BOOST_AUTO_TEST_CASE(NetworkModel_test_reaction_rule)
{
    Species sp1("A"), sp2("B"), sp3("C");

    ReactionRule rr1, rr2, rr3;
    rr1.add_reactant(sp1);
    rr1.add_reactant(sp2);
    rr1.add_product(sp3);
    rr2.add_reactant(sp3);
    rr2.add_product(sp1);
    rr2.add_product(sp2);
    rr3.add_reactant(sp1);
    rr3.add_product(sp2);

    NetworkModel model;
    model.add_reaction_rule(rr1);
    model.add_reaction_rule(rr2);
    BOOST_CHECK(model.has_reaction_rule(rr1));
    BOOST_CHECK(model.has_reaction_rule(rr2));
    BOOST_CHECK(!model.has_reaction_rule(rr3));
    model.add_reaction_rule(rr3);
    // BOOST_CHECK_THROW(model.add_reaction_rule(rr1), AlreadyExists); //XXX:
    model.remove_reaction_rule(rr1);
    BOOST_CHECK_THROW(model.remove_reaction_rule(rr1), NotFound);
    model.remove_reaction_rule(rr3);
    model.remove_reaction_rule(rr2);
}

BOOST_AUTO_TEST_CASE(NetworkModel_test_query_reaction_rules1)
{
    Species sp1("A"), sp2("B"), sp3("C");

    ReactionRule rr1, rr2, rr3, rr4;
    rr1.add_reactant(sp1);
    rr1.add_reactant(sp2);
    rr1.add_product(sp3);
    rr2.add_reactant(sp3);
    rr2.add_product(sp1);
    rr2.add_product(sp2);
    rr3.add_reactant(sp1);
    rr3.add_product(sp2);
    rr4.add_reactant(sp1);
    rr4.add_product(sp3);

    NetworkModel model;
    model.add_reaction_rule(rr1);
    model.add_reaction_rule(rr2);
    model.add_reaction_rule(rr3);
    model.add_reaction_rule(rr4);

    BOOST_CHECK_EQUAL(model.query_reaction_rules(sp1).size(), 2);
    BOOST_CHECK_EQUAL(model.query_reaction_rules(sp3).size(), 1);
    BOOST_CHECK_EQUAL(model.query_reaction_rules(sp2).size(), 0);
    BOOST_CHECK_EQUAL(model.query_reaction_rules(sp1, sp2).size(), 1);
    BOOST_CHECK((*(model.query_reaction_rules(sp1, sp2).begin())) == rr1);
}

BOOST_AUTO_TEST_CASE(NetworkModel_test_query_reaction_rules2)
{
    Species sp1("A"), sp2("B");

    ReactionRule rr1;
    rr1.add_reactant(sp1);
    rr1.add_reactant(sp2);

    NetworkModel model;
    model.add_reaction_rule(rr1);

    BOOST_CHECK_EQUAL(model.query_reaction_rules(sp1, sp2).size(), 1);
    BOOST_CHECK_EQUAL(model.query_reaction_rules(sp2, sp1).size(), 1);
}

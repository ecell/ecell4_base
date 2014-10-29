#define BOOST_TEST_MODULE "NetfreeModel_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/NetfreeModel.hpp>
#include <ecell4/core/extras.hpp>

using namespace ecell4;


BOOST_AUTO_TEST_CASE(NetfreeModel_test_constructor)
{
    NetfreeModel model;
}

BOOST_AUTO_TEST_CASE(NetfreeModel_test_species)
{
    Species sp1("A"), sp2("B");

    NetfreeModel model;
    model.add_species_attribute(sp1);
    model.add_species_attribute(sp2);
    BOOST_CHECK(model.has_species_attribute(sp1));
    BOOST_CHECK(model.has_species_attribute(sp2));
    BOOST_CHECK_THROW(model.add_species_attribute(sp1), AlreadyExists);
    model.remove_species_attribute(sp1);
    BOOST_CHECK_THROW(model.remove_species_attribute(sp1), NotFound);
}

BOOST_AUTO_TEST_CASE(NetfreeModel_test_reaction_rule)
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

    NetfreeModel model;
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

BOOST_AUTO_TEST_CASE(NetfreeModel_test_query_reaction_rules1)
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

    NetfreeModel model;
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

BOOST_AUTO_TEST_CASE(NetfreeModel_test_query_reaction_rules2)
{
    Species sp1("A"), sp2("B");

    ReactionRule rr1;
    rr1.add_reactant(sp1);
    rr1.add_reactant(sp2);

    NetfreeModel model;
    model.add_reaction_rule(rr1);

    BOOST_CHECK_EQUAL(model.query_reaction_rules(sp1, sp2).size(), 1);
    BOOST_CHECK_EQUAL(model.query_reaction_rules(sp2, sp1).size(), 1);
}

BOOST_AUTO_TEST_CASE(NetfreeModel_generation1)
{
    NetfreeModel nfm;
    nfm.add_reaction_rule(
        create_synthesis_reaction_rule(Species("X(p,q=a)"), 1.0));
    nfm.add_reaction_rule(
        create_unimolecular_reaction_rule(Species("X(q=a)"), Species("X(q=b)"), 1.0));
    nfm.add_reaction_rule(
        create_unbinding_reaction_rule(
            Species("X(p^1).X(p^1)"), Species("X(p)"), Species("X(p)"), 1.0));
    nfm.add_reaction_rule(
        create_binding_reaction_rule(
            Species("X(p)"), Species("X(p)"), Species("X(p^1).X(p^1)"), 1.0));

    std::vector<Species> seeds(1);
    seeds[0] = Species("X(p^1,q=a).X(p^1,q=a)");
    // seeds[1] = Species("X(p,q=a)");

    NetworkModel nwm(extras::generate_network_from_netfree_model(nfm, seeds, 10));

    for (NetworkModel::reaction_rule_container_type::const_iterator
        i(nwm.reaction_rules().begin()); i != nwm.reaction_rules().end(); ++i)
    {
        NetworkModel::reaction_rule_container_type::difference_type
            idx(std::distance(nwm.reaction_rules().begin(), i));
        std::cout << "[" << idx << "]: " << (*i).as_string() << std::endl;
    }

    BOOST_CHECK_EQUAL(nwm.reaction_rules().size(), 11);
}

BOOST_AUTO_TEST_CASE(NetfreeModel_generation2)
{
    NetfreeModel nfm;
    nfm.add_reaction_rule(
        create_synthesis_reaction_rule(Species("X(l,r)"), 1.0));
    nfm.add_reaction_rule(
        create_binding_reaction_rule(
            Species("X(r)"), Species("X(l)"), Species("X(r^1).X(l^1)"), 1.0));
    nfm.add_reaction_rule(
        create_unbinding_reaction_rule(
             Species("X(r^1).X(l^1)"),Species("X(r)"), Species("X(l)"), 1.0));

    std::vector<Species> seeds(0);
    utils::get_mapper_mf<Species, Integer>::type max_stoich;
    max_stoich[Species("X")] = 5;

    NetworkModel nwm(extras::generate_network_from_netfree_model(
        nfm, seeds, 2, max_stoich));

    for (NetworkModel::reaction_rule_container_type::const_iterator
        i(nwm.reaction_rules().begin()); i != nwm.reaction_rules().end(); ++i)
    {
        NetworkModel::reaction_rule_container_type::difference_type
            idx(std::distance(nwm.reaction_rules().begin(), i));
        std::cout << "[" << idx << "]: " << (*i).as_string() << std::endl;
    }

    BOOST_CHECK_EQUAL(nwm.reaction_rules().size(), 11);
}

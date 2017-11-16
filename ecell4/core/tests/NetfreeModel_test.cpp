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

    boost::shared_ptr<Model> nwm(nfm.expand(seeds, 10));

    // for (NetworkModel::reaction_rule_container_type::const_iterator
    //     i((*nwm).reaction_rules().begin()); i != (*nwm).reaction_rules().end(); ++i)
    // {
    //     NetworkModel::reaction_rule_container_type::difference_type
    //         idx(std::distance((*nwm).reaction_rules().begin(), i));
    //     std::cout << "[" << idx << "]: " << (*i).as_string() << std::endl;
    // }

    BOOST_CHECK_EQUAL((*nwm).reaction_rules().size(), 10);
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
    std::map<Species, Integer> max_stoich;
    max_stoich[Species("X")] = 5;

    boost::shared_ptr<Model> nwm(nfm.expand(seeds, 10, max_stoich));

    // for (NetworkModel::reaction_rule_container_type::const_iterator
    //     i((*nwm).reaction_rules().begin()); i != (*nwm).reaction_rules().end(); ++i)
    // {
    //     NetworkModel::reaction_rule_container_type::difference_type
    //         idx(std::distance((*nwm).reaction_rules().begin(), i));
    //     std::cout << "[" << idx << "]: " << (*i).as_string() << std::endl;
    // }

    BOOST_CHECK_EQUAL((*nwm).reaction_rules().size(), 13);
}

// BOOST_AUTO_TEST_CASE(NetfreeModel_query_reaction_rules3)
// {
// }

BOOST_AUTO_TEST_CASE(NetfreeModel_generation3)
{
    NetfreeModel m1;
    m1.add_reaction_rule(
        create_binding_reaction_rule(
            Species("A(r)"), Species("A(l)"), Species("A(r^1).A(l^1)"), 1.0));

    std::vector<ReactionRule> const retval1 = m1.query_reaction_rules(Species("A(l, r)"), Species("A(l, r)"));
    BOOST_CHECK_EQUAL(retval1.size(), 1);
    BOOST_CHECK_EQUAL(retval1[0].k(), 2.0);
    BOOST_CHECK_EQUAL(retval1[0].reactants().size(), 2);
    BOOST_CHECK_EQUAL(retval1[0].reactants()[0], Species("A(l,r)"));
    BOOST_CHECK_EQUAL(retval1[0].reactants()[1], Species("A(l,r)"));
    BOOST_CHECK_EQUAL(retval1[0].products().size(), 1);
    BOOST_CHECK_EQUAL(retval1[0].products()[0], Species("A(l,r^1).A(l^1,r)"));

    std::vector<Species> seeds1(1, Species("A(l, r)"));
    std::map<Species, Integer> max_stoich;
    max_stoich[Species("A")] = 4;
    boost::shared_ptr<Model> m2(m1.expand(seeds1, 100, max_stoich));
    std::vector<ReactionRule> const& reaction_rules1 = m2->reaction_rules();
    BOOST_CHECK_EQUAL(reaction_rules1.size(), 4);
    BOOST_CHECK_EQUAL(reaction_rules1[0].k(), 2.0);
    BOOST_CHECK_EQUAL(reaction_rules1[1].k(), 2.0);
    BOOST_CHECK_EQUAL(reaction_rules1[2].k(), 2.0);
    BOOST_CHECK_EQUAL(reaction_rules1[3].k(), 2.0);

    NetfreeModel m3;
    m3.add_reaction_rule(
        create_binding_reaction_rule(
            Species("A(r)"), Species("A(r)"), Species("A(r^1).A(r^1)"), 1.0));

    std::vector<ReactionRule> const retval2 = m3.query_reaction_rules(Species("A(r)"), Species("A(r)"));
    BOOST_CHECK_EQUAL(retval2.size(), 1);
    BOOST_CHECK_EQUAL(retval2[0].k(), 1.0);
    BOOST_CHECK_EQUAL(retval2[0].reactants().size(), 2);
    BOOST_CHECK_EQUAL(retval2[0].reactants()[0], Species("A(r)"));
    BOOST_CHECK_EQUAL(retval2[0].reactants()[1], Species("A(r)"));
    BOOST_CHECK_EQUAL(retval2[0].products().size(), 1);
    BOOST_CHECK_EQUAL(retval2[0].products()[0], Species("A(r^1).A(r^1)"));

    std::vector<ReactionRule> const retval3 = m3.query_reaction_rules(Species("A(r=u)"), Species("A(r=p)"));
    BOOST_CHECK_EQUAL(retval3.size(), 1);
    BOOST_CHECK_EQUAL(retval3[0].k(), 1.0);
    BOOST_CHECK_EQUAL(retval3[0].reactants().size(), 2);
    BOOST_CHECK_EQUAL(retval3[0].reactants()[0], Species("A(r=u)"));
    BOOST_CHECK_EQUAL(retval3[0].reactants()[1], Species("A(r=p)"));
    BOOST_CHECK_EQUAL(retval3[0].products().size(), 1);
    BOOST_CHECK_EQUAL(retval3[0].products()[0], Species("A(r=p^1).A(r=u^1)"));

    NetfreeModel m4;
    m4.add_reaction_rule(
        create_binding_reaction_rule(
            Species("_(b)"), Species("_(b)"), Species("_(b^1)._(b^1)"), 1.0));
    m4.add_reaction_rule(
        create_unbinding_reaction_rule(
            Species("_(b^1)._(b^1)"), Species("_(b)"), Species("_(b)"), 1.0));

    std::vector<ReactionRule> const retval4 = m4.query_reaction_rules(Species("A(b^1).A(b^1)"));
    BOOST_CHECK_EQUAL(retval4.size(), 1);
    BOOST_CHECK_EQUAL(retval4[0].k(), 1.0);
    BOOST_CHECK_EQUAL(retval4[0].reactants().size(), 1);
    BOOST_CHECK_EQUAL(retval4[0].reactants()[0], Species("A(b^1).A(b^1)"));
    BOOST_CHECK_EQUAL(retval4[0].products().size(), 2);
    BOOST_CHECK_EQUAL(retval4[0].products()[0], Species("A(b)"));
    BOOST_CHECK_EQUAL(retval4[0].products()[1], Species("A(b)"));

    std::vector<Species> seeds2(1, Species("A(b)"));
    boost::shared_ptr<Model> m5(m4.expand(seeds2));
    std::vector<ReactionRule> const& reaction_rules2 = m5->reaction_rules();
    BOOST_CHECK_EQUAL(reaction_rules2.size(), 2);

    std::vector<Species> seeds3(2);
    seeds3[0] = Species("A(b)");
    seeds3[1] = Species("B(b)");
    boost::shared_ptr<Model> m6(m4.expand(seeds3));
    std::vector<ReactionRule> const& reaction_rules3 = m6->reaction_rules();
    BOOST_CHECK_EQUAL(reaction_rules3.size(), 6);
    for (std::vector<ReactionRule>::const_iterator i(reaction_rules3.begin());
        i != reaction_rules3.end(); ++i)
    {
        BOOST_CHECK_EQUAL((*i).k(), 1.0);
    }
}

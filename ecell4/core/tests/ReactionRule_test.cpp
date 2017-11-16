#define BOOST_TEST_MODULE "ReactionRule_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/floating_point_comparison.hpp>

#include <ecell4/core/Species.hpp>
#include <ecell4/core/Context.hpp>
#include <ecell4/core/ReactionRule.hpp>

using namespace ecell4;


BOOST_AUTO_TEST_CASE(ReactionRule_test_constructor)
{
    ReactionRule rr;
}

BOOST_AUTO_TEST_CASE(ReactionRule_test_k)
{
    const Real epsrel(1e-6);
    const Real k(1.5);
    ReactionRule rr;
    rr.set_k(0);
    rr.set_k(k);
    BOOST_CHECK_CLOSE(rr.k(), k, k * epsrel);
    BOOST_CHECK_THROW(rr.set_k(-k), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(ReactionRule_test_reactants)
{
    ReactionRule rr;
    Species sp1("A"), sp2("B");
    rr.add_reactant(sp1);
    rr.add_reactant(sp2);
    rr.add_reactant(sp1);

    const ReactionRule::reactant_container_type& reactants(rr.reactants());
    BOOST_CHECK_EQUAL(reactants.size(), 3);
    BOOST_CHECK_EQUAL(std::count(reactants.begin(), reactants.end(), sp1), 2);
    BOOST_CHECK_EQUAL(std::count(reactants.begin(), reactants.end(), sp2), 1);
}

BOOST_AUTO_TEST_CASE(ReactionRule_test_products)
{
    ReactionRule rr;
    Species sp1("A"), sp2("B");
    rr.add_product(sp1);
    rr.add_product(sp2);
    rr.add_product(sp1);

    const ReactionRule::product_container_type& products(rr.products());
    BOOST_CHECK_EQUAL(products.size(), 3);
    BOOST_CHECK_EQUAL(std::count(products.begin(), products.end(), sp1), 2);
    BOOST_CHECK_EQUAL(std::count(products.begin(), products.end(), sp2), 1);
}

BOOST_AUTO_TEST_CASE(ReactionRule_test_compare)
{
    ReactionRule rr1, rr2, rr3, rr4, rr5;
    Species sp1("A"), sp2("B"), sp3("C"), sp4("D");

    rr1.add_reactant(sp1);
    rr1.add_reactant(sp2);
    rr1.add_product(sp3);
    rr1.set_k(1.5);

    rr2.add_reactant(sp1);
    rr2.add_reactant(sp3);
    rr2.add_product(sp3);
    rr2.set_k(1.5);

    rr3.add_reactant(sp1);
    rr3.add_reactant(sp2);
    rr3.add_product(sp4);
    rr3.set_k(1.5);

    rr4.add_reactant(sp2);
    rr4.add_reactant(sp1);
    rr4.add_product(sp3);
    rr4.set_k(5.0);

    rr5.add_reactant(sp1);
    rr5.add_reactant(sp2);
    rr5.add_product(sp3);
    rr5.set_k(5.0);

    BOOST_CHECK(rr1 == rr1);
    BOOST_CHECK(rr1 != rr2);
    BOOST_CHECK(rr1 != rr3);
    BOOST_CHECK(rr1 != rr4);
    BOOST_CHECK(rr1 == rr5);
}

void check_reaction_rule_generation(
    const ReactionRule& rr, const std::vector<Species>& reactants,
    const std::vector<ReactionRule>& ans,
    const std::string& filename, const unsigned int& lineno)
{
    const std::vector<ReactionRule> res = rr.generate(reactants);

    BOOST_CHECK_MESSAGE(res.size() == ans.size(), "" << filename << "(" << lineno << "): check res.size() == ans.size() failed [" << res.size() << " != " << ans.size() << "]");

    std::vector<std::vector<ReactionRule>::difference_type> done;
    done.reserve(res.size());
    for (std::vector<ReactionRule>::const_iterator i(res.begin());
        i != res.end(); ++i)
    {
        bool found = false;
        const ReactionRule formatted = format_reaction_rule_with_nosort(*i);
        for (std::vector<ReactionRule>::const_iterator j(ans.begin());
            j != ans.end(); ++j)
        {
            if (formatted == (*j) && formatted.k() == (*j).k())
            {
                std::vector<ReactionRule>::difference_type idx = std::distance(ans.begin(), j);
                if (std::find(done.begin(), done.end(), idx) == done.end())
                {
                    found = true;
                    done.push_back(idx);
                    break;
                }
            }
        }

        BOOST_CHECK_MESSAGE(found, "" << filename << "(" << lineno << "): A result [" << (*i).as_string() << "] is not in the list");
    }
}

#define ECELL4_TEST_REACTION_RULE_GENERATION( rr, reactants, ans ) \
    check_reaction_rule_generation(rr, reactants, ans, __FILE__, __LINE__)


BOOST_AUTO_TEST_CASE(ReactionRule_test_generate1)
{
    {
        const ReactionRule rr = create_synthesis_reaction_rule(Species("A"), 1.0);
        std::vector<Species> reactants;  // XXX: empty
        std::vector<ReactionRule> ans;
        ans.push_back(rr);

        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
}

BOOST_AUTO_TEST_CASE(ReactionRule_test_generate2)
{
    {
        const ReactionRule rr = create_binding_reaction_rule(
            Species("_1(b)"), Species("_1(b)"), Species("_1(b^1)._1(b^1)"), 1.0);

        ReactionRule::reactant_container_type reactants;
        reactants.push_back(Species("A(a^1,b).B(a^1,b)"));
        reactants.push_back(Species("B(a,b)"));

        std::vector<ReactionRule> ans;
        ans.push_back(format_reaction_rule_with_nosort(create_binding_reaction_rule(
            reactants[0], reactants[1], Species("A(a^1,b).B(a^1,b^3).B(a,b^3)"), 1.0)));

        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }

    {
        const ReactionRule rr = create_binding_reaction_rule(
            Species("A(b)"), Species("B(b)"), Species("A(b^1).B(b^1)"), 1.0);

        ReactionRule::reactant_container_type reactants;
        reactants.push_back(Species("A(a^1,b).A(a^1,b)"));
        reactants.push_back(Species("B(a^1,b).B(a^1,b^2).B(a^2,b)"));

        std::vector<ReactionRule> ans;
        const ReactionRule _ans1 = format_reaction_rule_with_nosort(create_binding_reaction_rule(
            reactants[0], reactants[1], Species("A(a^1,b^4).A(a^1,b).B(a^2,b^4).B(a^2,b^3).B(a^3,b)"), 1.0));
        ans.push_back(_ans1);
        ans.push_back(_ans1);
        const ReactionRule _ans2 = format_reaction_rule_with_nosort(create_binding_reaction_rule(
            reactants[0], reactants[1], Species("A(a^1,b^4).A(a^1,b).B(a^2,b).B(a^2,b^3).B(a^3,b^4)"), 1.0));
        ans.push_back(_ans2);
        ans.push_back(_ans2);

        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }

    std::vector<ReactionRule> retval;

    {
        const ReactionRule rr = create_unimolecular_reaction_rule(
            Species("A"), Species("B"), 1.0);
        ReactionRule::reactant_container_type reactants;
        reactants.push_back(Species("A"));
        std::vector<ReactionRule> ans;
        ans.push_back(rr);

        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }

    {
        const ReactionRule rr = create_unbinding_reaction_rule(
            Species("A(b^1).B(b^1)"), Species("A(b)"), Species("B(b)"), 1.0);
        ReactionRule::reactant_container_type reactants;
        reactants.push_back(
            Species("A(a^1,b^5).A(a^1,b^4).B(a^2,b).B(a^2,b^3).B(a^3,b^4).B(a,b^5)"));

        std::vector<ReactionRule> ans;
        ans.push_back(format_reaction_rule_with_nosort(create_unbinding_reaction_rule(
            reactants[0], Species("A(a^1,b).A(a^1,b^4).B(a^2,b).B(a^2,b^3).B(a^3,b^4)"), Species("B(a,b)"), 1.0)));
        ans.push_back(format_reaction_rule_with_nosort(create_unbinding_reaction_rule(
            reactants[0], Species("A(a^1,b^5).A(a^1,b).B(a,b^5)"), Species("B(a^2,b).B(a^2,b^3).B(a^3,b)"), 1.0)));

        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
}

BOOST_AUTO_TEST_CASE(ReactionRule_test_recursive_generation1)
{
    ReactionRule rr1;
    rr1.add_reactant(Species("X(r^1).X(l^1)"));
    rr1.add_product(Species("X(r)"));
    rr1.add_product(Species("X(l)"));
    rr1.set_k(1.0);
    const Species sp1("X(l,r^1).X(l^1,r^2).X(l^2,r^3).X(l^3,r^4).X(l^4,r)");

    ReactionRuleExpressionMatcher rrexp(rr1);
    BOOST_CHECK(rrexp.match(sp1));

    unsigned int i(0);
    do {
        ++i;
        std::vector<Species> products(rrexp.generate());
        // const ReactionRule tmp(rrexp.reactants(), products, rr1.k());
        // std::cerr << "GEN: " << tmp.as_string() << std::endl;
    } while (rrexp.next());

    BOOST_CHECK_EQUAL(i, 4);

    std::vector<ReactionRule> retval(rr1.generate(rrexp.reactants()));
    BOOST_CHECK_EQUAL(retval.size(), 4);
}

BOOST_AUTO_TEST_CASE(ReactionRule_test_generate3)
{
    {
        const ReactionRule rr = create_unbinding_reaction_rule(
            Species("A(b=u^1).A(b=u^1)"), Species("A(b=u)"), Species("A(b=u)"), 1.0);
        ReactionRule::reactant_container_type reactants(1, Species("A(b=u^1).A(b=u^1)"));
        std::vector<ReactionRule> ans;
        ans.push_back(format_reaction_rule_with_nosort(rr));
        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
    {
        const ReactionRule rr = create_unbinding_reaction_rule(
            Species("A(b=u^1).A(b=u^1)"), Species("A(b=u)"), Species("A(b=p)"), 1.0);
        ReactionRule::reactant_container_type reactants(1, Species("A(b=u^1).A(b=u^1)"));
        std::vector<ReactionRule> ans;
        ans.push_back(format_reaction_rule_with_nosort(rr));
        // ans.push_back(format_reaction_rule_with_nosort(rr));
        ans.push_back(format_reaction_rule_with_nosort(create_unbinding_reaction_rule(
            Species("A(b=u^1).A(b=u^1)"), Species("A(b=p)"), Species("A(b=u)"), 1.0)));
        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
    {
        const ReactionRule rr = create_unbinding_reaction_rule(
            Species("A(b^1,c^2).A(b^1,c^3).B(l^2,r^4).B(l^3,r^4)"),
            Species("A(b^1,c).A(b^1,c)"),
            Species("B(l,r^1).B(l,r^1)"), 1.0);
        ReactionRule::reactant_container_type reactants(1, Species("A(b^1,c^2).A(b^1,c^3).B(l^2,r^4).B(l^3,r^4)"));
        std::vector<ReactionRule> ans;
        ans.push_back(format_reaction_rule_with_nosort(create_unbinding_reaction_rule(
            reactants[0], Species("A(b^1,c).A(b^1,c)"), Species("B(l,r^1).B(l,r^1)"), 1.0)));
        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
    {
        ReactionRule rr = create_binding_reaction_rule(
            Species("A(b^1,c).A(b^1,c)"), Species("B(l,r^1).B(l,r^1)"),
            Species("A(b^1,c^2).A(b^1,c^3).B(l^2,r^4).B(l^3,r^4)"), 1.0);

        ReactionRule::reactant_container_type reactants(2);
        reactants[0] = Species("A(b^1,c).A(b^1,c)");
        reactants[1] = Species("B(l,r^1).B(l,r^1)");

        std::vector<ReactionRule> ans;
        ans.push_back(format_reaction_rule_with_nosort(rr));
        ans.push_back(format_reaction_rule_with_nosort(rr));

        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
    {
        ReactionRule rr = create_unbinding_reaction_rule(
            Species("_(b^1)._(b^1)"), Species("_(b)"), Species("_(b)"), 1.0);

        ReactionRule::reactant_container_type reactants(1, Species("A(b^1).A(b^1)"));
        std::vector<ReactionRule> ans;
        ans.push_back(format_reaction_rule_with_nosort(create_unbinding_reaction_rule(
            Species("A(b^1).A(b^1)"), Species("A(b)"), Species("A(b)"), 1.0)));

        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
    {
        ReactionRule rr = create_unbinding_reaction_rule(
            Species("_1(b^1)._1(b^1)"), Species("_1(b)"), Species("_1(b)"), 1.0);

        ReactionRule::reactant_container_type reactants(1, Species("A(b^1).B(b^1)"));
        std::vector<ReactionRule> ans;

        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
}

BOOST_AUTO_TEST_CASE(ReactionRule_test_generate4)
{
    {
        const ReactionRule rr = create_unimolecular_reaction_rule(
            Species("A(b^1).B(b^1)"), Species("A(b)"), 1.0);
        ReactionRule::reactant_container_type reactants(1, Species("A(b=u^1).B(b^1)"));
        std::vector<ReactionRule> ans;
        ans.push_back(format_reaction_rule_with_nosort(create_unimolecular_reaction_rule(
            reactants[0], Species("A(b=u)"), 1.0)));
        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
    {
        ReactionRule rr = create_degradation_reaction_rule(Species("B"), 1.0);
        rr.set_policy(ReactionRule::IMPLICIT);

        ReactionRule::reactant_container_type reactants(1, Species("A(b=u^1).B(b^1)"));
        std::vector<ReactionRule> ans;
        ans.push_back(format_reaction_rule_with_nosort(create_unimolecular_reaction_rule(
            reactants[0], Species("A(b=u)"), 1.0)));
        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
    {
        ReactionRule rr = create_degradation_reaction_rule(Species("B"), 1.0);
        rr.set_policy(ReactionRule::DESTROY);

        ReactionRule::reactant_container_type reactants(1, Species("A(b=u^1).B(b^1)"));
        std::vector<ReactionRule> ans;
        ans.push_back(format_reaction_rule_with_nosort(create_degradation_reaction_rule(reactants[0], 1.0)));
        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
    {
        const ReactionRule rr = create_unimolecular_reaction_rule(
            Species("A(b^1).A(b^1)"), Species("A(b)"), 1.0);
        ReactionRule::reactant_container_type reactants(1, Species("A(b^1).A(b^1)"));
        std::vector<ReactionRule> ans;
        const ReactionRule _ans = format_reaction_rule_with_nosort(create_unimolecular_reaction_rule(
            reactants[0], Species("A(b)"), 1.0));
        ans.push_back(_ans);
        ans.push_back(_ans);   //XXX: res should have two ReactionRules as shown here
        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
    {
        ReactionRule rr = create_degradation_reaction_rule(Species("A(b^1).A(b^1)"), 1.0);

        ReactionRule::reactant_container_type reactants(1, Species("A(b^1).A(b^1)"));
        std::vector<ReactionRule> ans;
        ans.push_back(format_reaction_rule_with_nosort(create_degradation_reaction_rule(reactants[0], 1.0)));
        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
    {
        const ReactionRule rr = create_unimolecular_reaction_rule(
            Species("A(b=_1^1).B(b^1)"), Species("A(b=_1^1).C(b=_1^1)"), 1.0);
        ReactionRule::reactant_container_type reactants(1, Species("A(b=u^1).B(b^1)"));
        std::vector<ReactionRule> ans;
        ans.push_back(format_reaction_rule_with_nosort(create_unimolecular_reaction_rule(
            reactants[0], Species("A(b=u^1).C(b=u^1)"), 1.0)));
        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
}

BOOST_AUTO_TEST_CASE(ReactionRule_test_generate5)
{
    {
        const ReactionRule rr = create_binding_reaction_rule(
            Species("X(r)"), Species("X(l)"), Species("X(r^1).X(l^1)"), 1.0);
        ReactionRule::reactant_container_type reactants(2, Species("X(l,r)"));
        std::vector<ReactionRule> ans;
        ans.push_back(format_reaction_rule_with_nosort(create_binding_reaction_rule(
            reactants[0], reactants[1], Species("X(l,r^1).X(l^1,r)"), 1.0)));
        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
    {
        const ReactionRule rr = create_binding_reaction_rule(
            Species("X(r)"), Species("X(l)"), Species("X(r^1).X(l^1)"), 1.0);
        ReactionRule::reactant_container_type reactants;
        reactants.push_back(Species("X(l,r^1).X(l^1,r)"));
        reactants.push_back(Species("X(l,r)"));
        std::vector<ReactionRule> ans;
        ans.push_back(format_reaction_rule_with_nosort(create_binding_reaction_rule(
            reactants[0], reactants[1], Species("X(l,r^1).X(l^1,r^2).X(l^2,r)"), 1.0)));
        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
    {
        const ReactionRule rr = create_unbinding_reaction_rule(
             Species("X(r^1).X(l^1)"),Species("X(r)"), Species("X(l)"), 1.0);
        ReactionRule::reactant_container_type reactants(1, Species("X(l,r^1).X(l^1,r)"));
        std::vector<ReactionRule> ans;
        ans.push_back(format_reaction_rule_with_nosort(create_unbinding_reaction_rule(
            reactants[0], Species("X(l,r)"), Species("X(l,r)"), 1.0)));
        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
    {
        const ReactionRule rr = create_unbinding_reaction_rule(
             Species("X(r^1).X(l^1)"),Species("X(r)"), Species("X(l)"), 1.0);
        ReactionRule::reactant_container_type reactants(1, Species("X(l,r^1).X(l^1,r^2).X(l^2,r)"));
        std::vector<ReactionRule> ans;
        ans.push_back(format_reaction_rule_with_nosort(create_unbinding_reaction_rule(
            reactants[0], Species("X(l,r^1).X(l^1,r)"), Species("X(l,r)"), 1.0)));
        ans.push_back(format_reaction_rule_with_nosort(create_unbinding_reaction_rule(
            reactants[0], Species("X(l,r)"), Species("X(l,r^1).X(l^1,r)"), 1.0)));
        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
}

BOOST_AUTO_TEST_CASE(ReactionRule_test_generate6)
{
    {
        ReactionRule rr = create_degradation_reaction_rule(Species("A"), 1.0);
        rr.set_policy(ReactionRule::IMPLICIT);
        ReactionRule::reactant_container_type reactants(1, Species("A(b)"));
        std::vector<ReactionRule> ans;
        const ReactionRule _ans = format_reaction_rule_with_nosort(
            create_degradation_reaction_rule(reactants[0], 1.0));
        ans.push_back(_ans);
        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
    {
        ReactionRule rr = create_degradation_reaction_rule(Species("A"), 1.0);
        rr.set_policy(ReactionRule::IMPLICIT);
        ReactionRule::reactant_container_type reactants(1, Species("A(b^1).B(b^1)"));
        std::vector<ReactionRule> ans;
        const ReactionRule _ans = format_reaction_rule_with_nosort(
            create_unimolecular_reaction_rule(reactants[0], Species("B(b)"), 1.0));
        ans.push_back(_ans);
        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
    {
        ReactionRule rr = create_degradation_reaction_rule(Species("A"), 1.0);
        rr.set_policy(ReactionRule::IMPLICIT);
        ReactionRule::reactant_container_type reactants(1, Species("A(b^1).A(b^1)"));
        std::vector<ReactionRule> ans;
        const ReactionRule _ans = format_reaction_rule_with_nosort(
            create_unimolecular_reaction_rule(reactants[0], Species("A(b)"), 1.0));
        ans.push_back(_ans);
        ans.push_back(_ans);
        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
    {
        ReactionRule rr = create_degradation_reaction_rule(Species("A(b^1).A(b^1)"), 1.0);
        ReactionRule::reactant_container_type reactants(1, Species("A(b^1).A(b^1)"));
        std::vector<ReactionRule> ans;
        const ReactionRule _ans = format_reaction_rule_with_nosort(
            create_degradation_reaction_rule(reactants[0], 1.0));
        ans.push_back(_ans);
        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
    {
        ReactionRule rr = create_degradation_reaction_rule(Species("A"), 1.0);
        rr.set_policy(ReactionRule::DESTROY);
        ReactionRule::reactant_container_type reactants(1, Species("A(b^1).A(b^1)"));
        std::vector<ReactionRule> ans;
        const ReactionRule _ans = format_reaction_rule_with_nosort(
            create_degradation_reaction_rule(reactants[0], 1.0));
        ans.push_back(_ans);
        ECELL4_TEST_REACTION_RULE_GENERATION(rr, reactants, ans);
    }
}

#undef ECELL4_TEST_REACTION_RULE_GENERATION

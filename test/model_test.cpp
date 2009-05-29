#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#define BOOST_TEST_MODULE "model"

#include <boost/test/included/unit_test.hpp>
#include <memory>
#include "Model.hpp"
#include "array_helper.hpp"
#include "SpeciesType.hpp"

BOOST_AUTO_TEST_CASE(basic)
{
    Model m;

    SpeciesType* s1 = m.new_species_type();
    (*s1)["name"] = "S";
    (*s1)["D"] = "1.5e-12";
    (*s1)["radius"] = "5e-9";

    SpeciesType* s2 = m.new_species_type();
    (*s2)["name"] = "P";
    (*s2)["D"] = "1e-12";
    (*s2)["radius"] = "7e-9";

    m.network_rules().add_reaction_rule(
        new_reaction_rule(s1, array_gen<SpeciesType*>(), .2));

    BOOST_CHECK_THROW(
        m.network_rules().add_reaction_rule(
            new_reaction_rule(s1, array_gen<SpeciesType*>(), .2)),
        already_exists);

    m.network_rules().add_reaction_rule(
        new_reaction_rule(s1, array_gen(s2), .2));

    BOOST_CHECK_THROW(
        m.network_rules().add_reaction_rule(
            new_reaction_rule(s1, array_gen(s2), .2)),
        already_exists);

    m.network_rules().add_reaction_rule(
        new_reaction_rule(s1, array_gen(s1, s2), .2));

    BOOST_CHECK_THROW(
        m.network_rules().add_reaction_rule(
            new_reaction_rule(s1, array_gen(s1, s2), .2)),
        already_exists);

    BOOST_CHECK_THROW(
        m.network_rules().add_reaction_rule(
            new_reaction_rule(s1, array_gen(s2, s1), .2)),
        already_exists);

    m.network_rules().add_reaction_rule(
        new_reaction_rule(s2, array_gen(s2, s1), .2));

    BOOST_CHECK_THROW(
        m.network_rules().add_reaction_rule(
            new_reaction_rule(s2, array_gen(s1, s2), .2)),
        already_exists);

    m.network_rules().add_reaction_rule(
        new_reaction_rule(s1, s2, array_gen(s2, s1), .2));

    BOOST_CHECK_THROW(
        m.network_rules().add_reaction_rule(
            new_reaction_rule(s1, s2, array_gen(s2, s1), .2)),
        already_exists);

}

BOOST_AUTO_TEST_CASE(query_reaction_rule)
{
    Model m;

    SpeciesType* s1 = m.new_species_type();
    (*s1)["name"] = "S";
    (*s1)["D"] = "1.5e-12";
    (*s1)["radius"] = "5e-9";

    SpeciesType* s2 = m.new_species_type();
    (*s2)["name"] = "P";
    (*s2)["D"] = "1e-12";
    (*s2)["radius"] = "7e-9";

    m.network_rules().add_reaction_rule(
        new_reaction_rule(s1, array_gen<SpeciesType*>(), .2));

    m.network_rules().add_reaction_rule(
        new_reaction_rule(s1, array_gen(s2), .2));

    m.network_rules().add_reaction_rule(
        new_reaction_rule(s2, array_gen<SpeciesType*>(), .2));

    {
        std::auto_ptr<NetworkRules::reaction_rule_generator> gen(
                m.network_rules().query_reaction_rule(s1));
        BOOST_CHECK(cue(*gen, new_reaction_rule(s1, array_gen<SpeciesType*>(), .2)));
    }

    {
        std::auto_ptr<NetworkRules::reaction_rule_generator> gen(
                m.network_rules().query_reaction_rule(s1));
        BOOST_CHECK(cue(*gen, new_reaction_rule(s1, array_gen(s2), .2)));
    }

    {
        std::auto_ptr<NetworkRules::reaction_rule_generator> gen(
                m.network_rules().query_reaction_rule(s1));
        BOOST_CHECK(!cue(*gen, new_reaction_rule(s1, array_gen(s1), .2)));
    }
}


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#define BOOST_TEST_MODULE "model"

#include <boost/test/included/unit_test.hpp>
#include <memory>
#include "Model.hpp"
#include "utils/array_helper.hpp"
#include "SpeciesType.hpp"

BOOST_AUTO_TEST_CASE(basic)
{
    Model m;

    boost::shared_ptr<SpeciesType> s1(new SpeciesType());
    (*s1)["name"] = "S";
    (*s1)["D"] = "1.5e-12";
    (*s1)["radius"] = "5e-9";

    boost::shared_ptr<SpeciesType> s2(new SpeciesType());
    (*s2)["name"] = "P";
    (*s2)["D"] = "1e-12";
    (*s2)["radius"] = "7e-9";

    m.add_species_type(s1);
    m.add_species_type(s2);

    m.network_rules().add_reaction_rule(
        new_reaction_rule(s1->id(), array_gen<SpeciesTypeID>(), .2));

    BOOST_CHECK_THROW(
        m.network_rules().add_reaction_rule(
            new_reaction_rule(s1->id(), array_gen<SpeciesTypeID>(), .2)),
        already_exists);

    m.network_rules().add_reaction_rule(
        new_reaction_rule(s1->id(), array_gen(s2->id()), .2));

    BOOST_CHECK_THROW(
        m.network_rules().add_reaction_rule(
            new_reaction_rule(s1->id(), array_gen(s2->id()), .2)),
        already_exists);

    m.network_rules().add_reaction_rule(
        new_reaction_rule(s1->id(), array_gen(s1->id(), s2->id()), .2));

    BOOST_CHECK_THROW(
        m.network_rules().add_reaction_rule(
            new_reaction_rule(s1->id(), array_gen(s1->id(), s2->id()), .2)),
        already_exists);

    BOOST_CHECK_THROW(
        m.network_rules().add_reaction_rule(
            new_reaction_rule(s1->id(), array_gen(s2->id(), s1->id()), .2)),
        already_exists);

    m.network_rules().add_reaction_rule(
        new_reaction_rule(s2->id(), array_gen(s2->id(), s1->id()), .2));

    BOOST_CHECK_THROW(
        m.network_rules().add_reaction_rule(
            new_reaction_rule(s2->id(), array_gen(s1->id(), s2->id()), .2)),
        already_exists);

    m.network_rules().add_reaction_rule(
        new_reaction_rule(s1->id(), s2->id(), array_gen(s2->id(), s1->id()), .2));

    BOOST_CHECK_THROW(
        m.network_rules().add_reaction_rule(
            new_reaction_rule(s1->id(), s2->id(), array_gen(s2->id(), s1->id()), .2)),
        already_exists);

}

BOOST_AUTO_TEST_CASE(query_reaction_rule)
{
    Model m;

    boost::shared_ptr<SpeciesType> s1(new SpeciesType());
    (*s1)["name"] = "S";
    (*s1)["D"] = "1.5e-12";
    (*s1)["radius"] = "5e-9";

    m.add_species_type(s1);

    boost::shared_ptr<SpeciesType> s2(new SpeciesType());
    (*s2)["name"] = "P";
    (*s2)["D"] = "1e-12";
    (*s2)["radius"] = "7e-9";

    m.add_species_type(s2);

    m.network_rules().add_reaction_rule(
        new_reaction_rule(s1->id(), array_gen<SpeciesTypeID>(), .2));

    m.network_rules().add_reaction_rule(
        new_reaction_rule(s1->id(), array_gen(s2->id()), .2));

    m.network_rules().add_reaction_rule(
        new_reaction_rule(s2->id(), array_gen<SpeciesTypeID>(), .2));

    {
        std::auto_ptr<NetworkRules::reaction_rule_generator> gen(
                m.network_rules().query_reaction_rule(s1->id()));
        BOOST_CHECK(cue(*gen, new_reaction_rule(s1->id(), array_gen<SpeciesTypeID>(), .2)));
    }

    {
        std::auto_ptr<NetworkRules::reaction_rule_generator> gen(
                m.network_rules().query_reaction_rule(s1->id()));
        BOOST_CHECK(cue(*gen, new_reaction_rule(s1->id(), array_gen(s2->id()), .2)));
    }

    {
        std::auto_ptr<NetworkRules::reaction_rule_generator> gen(
                m.network_rules().query_reaction_rule(s1->id()));
        BOOST_CHECK(!cue(*gen, new_reaction_rule(s1->id(), array_gen(s1->id()), .2)));
    }
}


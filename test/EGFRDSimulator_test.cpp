#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#define BOOST_TEST_MODULE "EGFRDSimulator"

#include <boost/test/included/unit_test.hpp>
#include "Model.hpp"
#include "EGFRDSimulator.hpp"
#include "NetworkRulesWrapper.hpp"
#include "ReactionRuleInfo.hpp"

typedef World<CyclicWorldTraits<Real, Real> > world_type;
typedef EGFRDSimulatorTraitsBase<world_type> simulator_traits_type;
typedef EGFRDSimulator<simulator_traits_type> simulator_type;
typedef simulator_traits_type::network_rules_type network_rules_type;

struct ReactionRecorderMock: public simulator_traits_type::reaction_recorder_type
{
    virtual ~ReactionRecorderMock() {}

    virtual void operator()(reaction_record_type const&) {}
};

BOOST_AUTO_TEST_CASE(test)
{

    boost::shared_ptr<world_type> w(new world_type(1., 10));
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

    boost::shared_ptr<network_rules_type> nrw(new network_rules_type(m.network_rules()));
    world_type::traits_type::rng_type rng;
    ReactionRecorderMock rr;
    simulator_type s(w, nrw, rng, rr, 3);

    s.step();
}

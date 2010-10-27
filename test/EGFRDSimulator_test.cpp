#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#define BOOST_TEST_MODULE "EGFRDSimulator"

#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include "ParticleModel.hpp"
#include "EGFRDSimulator.hpp"
#include "Logger.hpp"
#include "NullLogger.hpp"
#include "EGFRDSimulatorFactory.hpp"
#include "NetworkRulesWrapper.hpp"
#include "ReactionRuleInfo.hpp"

typedef World<CyclicWorldTraits<Real, Real> > world_type;
typedef EGFRDSimulatorTraitsBase<world_type> simulator_traits_type;
typedef EGFRDSimulator<simulator_traits_type> simulator_type;
typedef simulator_traits_type::network_rules_type network_rules_type;

static Real const N_A(6.0221367e23);

template<typename Tworld_, typename Trng_, typename Tpid_list_>
void inject_particles(Tworld_& world, Trng_& rng, Tpid_list_& pid_list, typename Tworld_::species_id_type const& sid, int n)
{
    typedef typename Tworld_::particle_id_pair_and_distance_list particle_id_pair_list;
    typedef typename Tworld_::particle_id_pair particle_id_pair;
    typedef typename Tworld_::length_type length_type;
    typedef typename Tworld_::particle_shape_type particle_shape_type;
    typedef typename Tworld_::position_type position_type;
    typedef typename Tworld_::species_type species_type;

    species_type const& s(world.get_species(sid));
 
    for (int i = 0; i < n; ++i)
    {
        particle_shape_type p(position_type(), s.radius());

        for (;;)
        {
            p.position() = position_type(
                rng.uniform(0, world.world_size()),
                rng.uniform(0, world.world_size()),
                rng.uniform(0, world.world_size()));
            if (boost::scoped_ptr<particle_id_pair_list>(
                world.check_overlap(p)) == 0)
            {
                break;
            }
            std::cerr << i << "th particle rejected" << std::endl;
        }
        
        particle_id_pair i(world.new_particle(sid, p.position()));
        pid_list.push_back(i.first);
    }
}

BOOST_AUTO_TEST_CASE(test)
{
    typedef world_type::traits_type::particle_id_type particle_id_type;
    typedef world_type::traits_type::length_type length_type;

    int const n(300);

    LoggerFactory::register_logger_factory(
        "EGFRDSimulator",
        boost::shared_ptr<LoggerFactory>(new NullLoggerFactory()));

    ParticleModel m;

    m["size"] = "1.";
    m["matrix_size"] = boost::lexical_cast<std::string>(
        (int)std::pow(n * 6, 1. / 3.));
    std::cout << "matrix_size=" << m["matrix_size"] << std::endl;

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
        new_reaction_rule(s1->id(), s1->id(),
            array_gen<SpeciesTypeID>(s2->id()), 1e7 / N_A));

    m.network_rules().add_reaction_rule(
        new_reaction_rule(s2->id(),
            array_gen<SpeciesTypeID>(s1->id(), s1->id()), 1e3));

    boost::shared_ptr<network_rules_type> nrw(
        new network_rules_type(m.network_rules()));

    world_type::traits_type::rng_type rng;

    std::vector<particle_id_type> S_particles, P_particles;
    EGFRDSimulatorFactory<simulator_traits_type> factory(rng);
    boost::scoped_ptr<simulator_type> s(factory(m));

    inject_particles(*s->world(), rng, S_particles, s1->id(), n / 2);
    inject_particles(*s->world(), rng, P_particles, s2->id(), n / 2);

    for (int i = 10000; --i >= 0;)
        s->step();

    BOOST_TEST_MESSAGE("spherical single: " << s->num_domains_per_type(simulator_type::SPHERICAL_SINGLE));
    BOOST_TEST_MESSAGE("cylindrical single: " << s->num_domains_per_type(simulator_type::CYLINDRICAL_SINGLE));
    BOOST_TEST_MESSAGE("spherical pair: " << s->num_domains_per_type(simulator_type::SPHERICAL_PAIR));
    BOOST_TEST_MESSAGE("cylindrical pair: " << s->num_domains_per_type(simulator_type::CYLINDRICAL_PAIR));
    BOOST_TEST_MESSAGE("multi: " << s->num_domains_per_type(simulator_type::MULTI));
}

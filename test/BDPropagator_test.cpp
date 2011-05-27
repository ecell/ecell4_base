#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#define BOOST_TEST_MODULE "BDPropagator"

#include <boost/test/included/unit_test.hpp>

#include "utils/range.hpp"
#include "utils/pair.hpp"
#include "MatrixSpace.hpp"
#include "utils.hpp"
#include "Sphere.hpp"
#include "Cylinder.hpp"
#include "Box.hpp"
#include "Surface.hpp"
#include "Region.hpp"
#include "Plane.hpp"
#include "BDPropagator.hpp"
#include "NetworkRules.hpp"
#include "Transaction.hpp"
#include "World.hpp"
#include "GSLRandomNumberGenerator.hpp"
#include "BasicNetworkRulesImpl.hpp"
#include "NetworkRulesWrapper.hpp"
#include "ReactionRuleInfo.hpp"
#include "ParticleSimulator.hpp"
#include "ReactionRecorder.hpp"
#include "VolumeClearer.hpp"

struct Traits: ParticleSimulatorTraitsBase<World<CyclicWorldTraits<Real, Real> > >
{};

typedef ParticleSimulator<Traits> _ParticleSimulator;

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

template<typename reaction_record_type>
struct reaction_recorder: ReactionRecorder<reaction_record_type>
{
    reaction_recorder() {}
    virtual ~reaction_recorder() {}
    virtual void operator()(reaction_record_type const& rec) {}
};

template<typename particle_shape_type, typename particle_id_type>
struct volume_clearer: VolumeClearer<particle_shape_type, particle_id_type>
{
    volume_clearer() {}
    virtual ~volume_clearer() {}

    virtual bool operator()(particle_shape_type const& shape, particle_id_type const& ignore)
    {
        return true;
    }

    virtual bool operator()(particle_shape_type const& shape, particle_id_type const& ignore0, particle_id_type const& ignore1)
    {
        return true;
    }
};

BOOST_AUTO_TEST_CASE(instantiation)
{
    Traits::world_type::traits_type::rng_type rng;
    BasicNetworkRulesImpl nr;
    Traits::network_rules_type nrw(nr);
    Traits::world_type w;
    boost::scoped_ptr<Traits::world_type::transaction_type> tx(
            w.create_transaction());
    reaction_recorder<Traits::reaction_record_type> rr;
    volume_clearer<Traits::world_type::particle_shape_type, Traits::world_type::particle_id_type> vc;
    BDPropagator<Traits> bdp(
        *tx, nrw, rng, 0.01, 100, &rr, &vc, 
        make_select_first_range(w.get_particles_range()));
}

BOOST_AUTO_TEST_CASE(basic)
{
    typedef Traits::world_type::species_id_type species_id;
    typedef Traits::world_type::particle_id_type particle_id;
    typedef Traits::world_type::species_type species;
    typedef Traits::world_type::particle_id_pair_generator particle_id_pair_generator;
    typedef Traits::network_rules_type::reaction_rule_type reaction_rule_type;
    typedef Traits::world_type::position_type position_type;
    SerialIDGenerator<species_id> sidgen;
    Traits::world_type::traits_type::rng_type rng;
    BasicNetworkRulesImpl nr;
    Traits::network_rules_type nrw(nr);
    Traits::world_type w(1e-5, 10);

    species S0(sidgen(), 2e-11, 5e-8, "default");
    species S1(sidgen(), 0, 1e-8, "default");
    species S2(sidgen(), 2e-11, 5e-8, "default");
    w.add_species(S0);
    w.add_species(S1);
    w.add_species(S2);

    boost::shared_ptr<Traits::world_type::structure_type> default_surface(
        new _ParticleSimulator::cuboidal_region_type("default",
            _ParticleSimulator::box_type(position_type(1e-5 / 2, 1e-5 / 2, 1e-5 / 2),
                             array_gen(1e-5, 1e-5, 1e-5))));
    w.add_structure(default_surface);

    nr.add_reaction_rule(new_reaction_rule(S0.id(), S1.id(), array_gen(S2.id()), 1e-9));

    std::vector<particle_id> S0_particles, S1_particles;
    inject_particles(w, rng, S0_particles, S0.id(), 500);
    inject_particles(w, rng, S1_particles, S1.id(), 500);
    int num_of_moving_particles = 500;
    int num_of_immobile_particles = 500;

    reaction_recorder<Traits::reaction_record_type> rr;
    volume_clearer<Traits::world_type::particle_shape_type, Traits::world_type::particle_id_type> vc;

    for (int i = 1000; --i >= 0; ) {
        boost::scoped_ptr<Traits::world_type::transaction_type> tx(w.create_transaction());
        BDPropagator<Traits> prpg(*tx, nrw, rng, 5e-11, 100, &rr, &vc, make_select_first_range(w.get_particles_range()));
        while (prpg());
        boost::scoped_ptr<particle_id_pair_generator> added_particles(tx->get_added_particles());
        boost::scoped_ptr<particle_id_pair_generator> removed_particles(tx->get_removed_particles());
        boost::scoped_ptr<particle_id_pair_generator> modified_particles(tx->get_modified_particles());
        BOOST_CHECK(count(*removed_particles) == 0 || count(*removed_particles) == 2);
        BOOST_CHECK_EQUAL(num_of_moving_particles, count(*removed_particles) / 2 + count(*modified_particles) + prpg.get_rejected_move_count());
        num_of_immobile_particles -= count(*removed_particles) / 2;
        BOOST_CHECK_EQUAL(w.num_particles(), num_of_moving_particles + num_of_immobile_particles);
    }
}


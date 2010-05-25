#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/timer.hpp>
#include <algorithm>

#include "utils/range.hpp"
#include "MatrixSpace.hpp"
#include "utils.hpp"
#include "Sphere.hpp"
#include "Cylinder.hpp"
#include "Box.hpp"
#include "Surface.hpp"
#include "Region.hpp"
#include "Plane.hpp"
#include "BDSimulator.hpp"
#include "NetworkRules.hpp"
#include "Transaction.hpp"
#include "World.hpp"
#include "GSLRandomNumberGenerator.hpp"
#include "BasicNetworkRulesImpl.hpp"
#include "NetworkRulesWrapper.hpp"
#include "ReactionRuleInfo.hpp"
#include "ParticleSimulator.hpp"

struct Traits: ParticleSimulatorTraitsBase<World<CyclicWorldTraits<Real, Real> > >
{};

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

void do_benchmark(Real volume, std::size_t n, Traits::time_type t, Real dt_factor)
{
    typedef Traits::world_type::species_id_type species_id;
    typedef Traits::world_type::particle_id_type particle_id;
    typedef Traits::world_type::species_type species;
    typedef Traits::world_type::length_type length_type;
    typedef Traits::world_type::particle_id_pair_generator particle_id_pair_generator;
    typedef Traits::network_rules_type::reaction_rule_type reaction_rule_type;
    typedef Traits::world_type::position_type position_type;
    SerialIDGenerator<species_id> sidgen;
    Traits::world_type::traits_type::rng_type rng;
    BasicNetworkRulesImpl nr;
    Traits::network_rules_type nrw(nr);
    length_type const world_size(std::pow(volume, 1. / 3.));
    std::size_t const matrix_size(
        std::max(static_cast<std::size_t>(3u),
                 static_cast<std::size_t>(std::pow(3. * n, 1. / 3.))));
    Traits::world_type w(world_size, matrix_size);

    species A(sidgen(), 1e-12, 2.5e-9, "default");
    w.add_species(A);

    boost::shared_ptr<Traits::world_type::structure_type> default_surface(
        new Traits::cuboidal_region_type("default",
            Traits::box_type(
                position_type(world_size / 2, world_size / 2, world_size / 2),
                 array_gen(world_size, world_size, world_size))));
    w.add_structure(default_surface);

    std::vector<particle_id> A_particles;
    inject_particles(w, rng, A_particles, A.id(), n);

    std::cout << "T: " << t << std::endl;
    std::cout << "V: " << volume << std::endl;
    std::cout << "N: " << n << std::endl;
    std::cout << "world size: " << world_size << std::endl;
    std::cout << "matrix size: " << matrix_size << std::endl;

    {
        std::cout << "stir" << std::endl;
        BDSimulator<Traits> s(w, nrw, rng, dt_factor);
        std::cout << "dt: " << s.dt() << std::endl;
        Traits::time_type const stir_time(t * .1);
        while (s.step(stir_time));
    }

    {
        std::cout << "run" << std::endl;
        BDSimulator<Traits> s(w, nrw, rng, dt_factor);
        boost::timer timer;
        while (s.step(t));
        std::cout << "t: " << s.t() << "=" << t << std::endl;
        std::cout << "dt: " << s.dt() << std::endl;
        std::cout << "steps (total): " << s.num_steps() << std::endl;
        std::cout << "elapsed: " << timer.elapsed() << std::endl;
        std::cout << "steps per second: "
                  << (static_cast<double>(s.num_steps()) / timer.elapsed())
                  << std::endl;
        std::cout << "steps/N: "
                  << (static_cast<double>(s.num_steps()) / n)
                  << std::endl;
    }
}

int main()
{
    do_benchmark(1e-12, 1e5, 1e-10, 1e-6);
    return 0;
}

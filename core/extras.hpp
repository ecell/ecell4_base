#ifndef __ECELL4_EXTRAS_HPP
#define __ECELL4_EXTRAS_HPP

#include "types.hpp"
#include "Position3.hpp"
#include "Species.hpp"
#include "Particle.hpp"


namespace ecell4
{

namespace extras
{

template<typename Tworld_, typename Trng_>
void throw_in_particles(
    Tworld_& world, const Species& sp, const Integer& N, Trng_& rng)
{
    typedef typename Tworld_::particle_info_type particle_info_type;

    if (N < 0)
    {
        throw std::invalid_argument("the number of particles must be positive.");
    }

    const Position3 edge_lengths(world.edge_lengths());
    const particle_info_type info(world.get_particle_info(sp));

    for (int i(0); i < N; ++i)
    {
        while (true)
        {
            const Position3 pos(
                rng.uniform(0.0, edge_lengths[0]),
                rng.uniform(0.0, edge_lengths[1]),
                rng.uniform(0.0, edge_lengths[2]));
            if (world.list_particles_within_radius(pos, info.radius).size()
                == 0)
            {
                world.new_particle(Particle(sp, pos, info.radius, info.D));
                break;
            }
        }
    }
}

} // extras

} // ecell4

#endif // __ECELL4_EXTRAS_HPP

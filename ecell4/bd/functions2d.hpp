#ifndef ECELL4_BD_FUNCTIONS_2D_HPP
#define ECELL4_BD_FUNCTIONS_2D_HPP

#include <ecell4/core/config.h>
#include <ecell4/core/types.hpp>
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include "rotate_vector.hpp"

namespace ecell4
{

namespace bd
{

inline Real3
random_circular_uniform(RandomNumberGenerator& rng, const Real& r)
{
    const Real theta = rng.uniform(0., 2 * M_PI);
    return Real3(r * std::cos(theta), r * std::sin(theta), 0.);
}

inline Real3
random_circular_uniform(RandomNumberGenerator& rng,
                        const Real& r, const Real3& normal)
{
    const Real theta = rng.uniform(0., 2 * M_PI);
    const Real l = r / std::sqrt(normal[0] * normal[0] + normal[1] * normal[1]);
    const Real3 retval(l * normal[1], -l * normal[0], 0.);
    return rotate(theta, normal, retval);
}

inline Real3
random_displacement_2d(RandomNumberGenerator& rng, const Real& t, const Real& D)
{
    const Real sigma(std::sqrt(2 * D * t));
    return Real3(
        rng.gaussian(sigma), rng.gaussian(sigma), 0);
}

inline Real3
random_displacement_2d(RandomNumberGenerator& rng, const Real& t, const Real& D,
                       const Real3& normal)
{
    const Real sigma(std::sqrt(4 * D * t));
    const Real len(rng.gaussian(sigma));
    return random_circular_uniform(rng, len, normal);
}

inline Real3 random_ipv_2d(RandomNumberGenerator& rng,
        const Real r, const Real D, const Real reaction_length)
{
    const Real rl    = r + reaction_length;
    const Real r_sq  = r * r;
    const Real rl_sq = rl * rl;

    const Real rnd = rng.uniform(0., 1.);
    const Real ipvl = std::sqrt(r_sq + rnd * (rl_sq - r_sq));

    return random_circular_uniform(rng, ipvl);
}

inline Real3 random_ipv_2d(RandomNumberGenerator& rng,
        const Real r, const Real D, const Real reaction_length,
        const Real3& normal)
{
    const Real rl    = r + reaction_length;
    const Real r_sq  = r * r;
    const Real rl_sq = rl * rl;

    const Real rnd = rng.uniform(0., 1.);
    const Real ipvl = std::sqrt(r_sq + rnd * (rl_sq - r_sq));

    return random_circular_uniform(rng, ipvl, normal);
}


} // bd

} // ecell4

#endif /* ECELL4_BD_FUNCTIONS_2D_HPP */

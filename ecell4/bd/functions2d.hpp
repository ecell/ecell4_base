#ifndef ECELL4_BD_FUNCTIONS_2D_HPP
#define ECELL4_BD_FUNCTIONS_2D_HPP

#include <ecell4/core/config.h>
#include <ecell4/core/types.hpp>
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/geometry.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

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
    const Real3 rnd = random_circular_uniform(rng, r);
    const Real3 unitz(0, 0, 1);
    const Real tilt = angle(unitz, normal);

    if(std::abs(tilt - M_PI) < 1e-12)      return rnd;
    else if(std::abs(tilt + M_PI) < 1e-12) return rnd * (-1.0);

    const Real3 cross = cross_product(unitz, normal);
    const Real3 retval = rotate(tilt, cross, rnd);
    return retval;
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

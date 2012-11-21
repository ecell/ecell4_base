#ifndef __FUNCTIONS_3D_HPP
#define __FUNCTIONS_3D_HPP

#include <ecell4/core/types.hpp>
#include <ecell4/core/Position3.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

namespace ecell4
{

namespace bd
{

Real I_bd_3d(Real const& sigma, Real const& t, Real const& D);

Position3 random_spherical_uniform(RandomNumberGenerator& rng, Real const& r);
Position3 random_displacement_3d(
    RandomNumberGenerator& rng, Real const& t, Real const& D);
Position3 random_ipv_3d(
    RandomNumberGenerator& rng, Real const& sigma, Real const& t, Real const& D);

} // bd

} // ecell4

#endif /* __FUNCTIONS_3D_HPP */

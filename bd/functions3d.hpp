#ifndef __FUNCTIONS_3D_HPP
#define __FUNCTIONS_3D_HPP

#include <ecell4/core/types.hpp>
#include <ecell4/core/Position3.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

namespace ecell4
{

namespace bd
{

Position3 random_unit_vector_3d(RandomNumberGenerator& rng);
Position3 random_displacement_3d(
    RandomNumberGenerator& rng, Real const& t, Real const& D);
Real drawR_gbd(
    RandomNumberGenerator& rng, Real const sigma, Real const& t, Real const& D);
Real I_bd_3d(Real const& sigma, Real const& t, Real const& D);

} // bd

} // ecell4

#endif /* __FUNCTIONS_3D_HPP */

#ifndef __FUNCTIONS_3D_HPP
#define __FUNCTIONS_3D_HPP

#include <ecell4/core/types.hpp>
#include <ecell4/core/Position3.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

namespace ecell4
{

namespace bd
{

/**
 * $\int_0^\infty r^2dr\,g\left(r,\Delta t\right),$
 * where $g\left(r,\Delta t\right)$ is a probability that a pair, which is
 * initially separated by a length $r$, overlaps after the time $\Delta t$:
 * $g\left(r,\Delta t\right)\equiv\int_0^\sigma r'^2dr'\int_0^\pi\sin\theta
 * d\theta\int_0^{2\pi}d\phi\,p\left({\bf r'},t+\Delta t;{\bf r},t\right).$
 * see Eqs. (20-21) in (Morelli & ten Wolde, J. Chem. Phys., 2008).
 * @param sigma a radius of the excluded volume, $\sigma$.
 * @param t a step interval, $\Delta t$.
 * @param D a diffusion coefficient, $D$.
 */
Real Igbd_3d(Real const& sigma, Real const& t, Real const& D);

/**
 * $\int_0^R r^2dr\,g\left(r,\Delta t\right).$
 * see Eqs. (20-21) in (Morelli & ten Wolde, J. Chem. Phys., 2008).
 * @param an upper limit of the integration, $R$.
 * @param sigma a radius of the excluded volume, $\sigma$.
 * @param t a step interval, $\Delta t$.
 * @param D a diffusion coefficient, $D$.
 */
Real Igbd_r_3d(Real r, Real sigma, Real t, Real D);

Position3 random_spherical_uniform(RandomNumberGenerator& rng, Real const& r);
Position3 random_displacement_3d(
    RandomNumberGenerator& rng, Real const& t, Real const& D);
Position3 random_ipv_3d(
    RandomNumberGenerator& rng, Real const& sigma, Real const& t, Real const& D);

} // bd

} // ecell4

#endif /* __FUNCTIONS_3D_HPP */

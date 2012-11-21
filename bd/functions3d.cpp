#include "functions3d.hpp"


namespace ecell4
{

namespace bd
{

Position3 random_unit_vector_3d(RandomNumberGenerator& rng)
{
    return Position3(1, 0, 0);
}

Position3 random_displacement_3d(
    RandomNumberGenerator& rng, Real const& t, Real const& D)
{
    Real const sigma(std::sqrt(2 * D * t));
    return Position3(
        rng.gaussian(0, sigma), rng.gaussian(0, sigma), rng.gaussian(0, sigma));
}

Real I_bd_3d(Real const& sigma, Real const& t, Real const& D)
{
    const Real sqrtPi(std::sqrt(M_PI));

    const Real Dt(D * t);
    const Real Dt2(Dt + Dt);
    const Real sqrtDt(std::sqrt(Dt));
    const Real sigmasq(sigma * sigma);

    const Real term1(1 / (3 * sqrtPi));
    const Real term2(sigmasq - Dt2);
    const Real term3(Dt2 - 3 * sigmasq);
    const Real term4(sqrtPi * sigmasq * sigma * erfc(sigma / sqrtDt));

    const Real result(
        term1 * (-sqrtDt * (term2 * std::exp(-sigmasq / Dt) + term3) + term4));
    return result;
}

Real drawR_gbd(
    RandomNumberGenerator& rng, Real const sigma, Real const& t, Real const& D)
{
    return 0;
}

} // bd

} // ecell4

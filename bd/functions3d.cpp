#include <gsl/gsl_roots.h>

#include "functions3d.hpp"


namespace ecell4
{

namespace bd
{

Position3 random_spherical_uniform(
    RandomNumberGenerator& rng, Real const& r)
{
    Real a(0), b(0), r2(1);
    while (r2 > 0.25)
    {
        a = rng.uniform(0, 1) - 0.5;
        b = rng.uniform(0, 1) - 0.5;
        r2 = a * a + b * b;
    }

    const Real scale(8 * r * std::sqrt(0.25 - r2));
    return Position3(a * scale, b * scale, r * (8 * r2 - 1));
}

Position3 random_displacement_3d(
    RandomNumberGenerator& rng, Real const& t, Real const& D)
{
    Real const sigma(std::sqrt(2 * D * t));
    return Position3(
        rng.gaussian(0, sigma), rng.gaussian(0, sigma), rng.gaussian(0, sigma));
}

Real Igbd_3d(Real const& sigma, Real const& t, Real const& D)
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

Real Igbd_r_3d(Real r, Real sigma, Real t, Real D)
{
    const Real sqrtPi(std::sqrt(M_PI));

    const Real Dt(D * t);
    const Real Dt2(Dt + Dt);
    const Real Dt4(Dt2 + Dt2);
    const Real sqrtDt(std::sqrt(Dt));
    // const Real sqrtDt4(std::sqrt(Dt4));
    const Real sqrtDt4(2 * sqrtDt);
    const Real sigmasq(sigma * sigma);

    const Real sigmacb(sigmasq * sigma);
    const Real rcb(gsl_pow_3(r));

    const Real rsigma(r * sigma);

    const Real rps_sq(gsl_pow_2(r + sigma)), rms_sq(gsl_pow_2(r - sigma));

    const Real term1(-2 * sqrtDt / sqrtPi);
    const Real term2(std::exp(-sigmasq / Dt) * (sigmasq - Dt2));
    const Real term3(-std::exp(-rps_sq / Dt4) * (rms_sq + rsigma - Dt2));
    const Real term4(std::exp(-rms_sq / Dt4) * (rps_sq - rsigma - Dt2));
    const Real term5(-sigmasq * 3 + Dt2);

    const Real term6((sigmacb - rcb) * erf((r - sigma) / sqrtDt4));
    const Real term7(-(sigmacb + sigmacb) * erf(sigma / sqrtDt));
    const Real term8((sigmacb + rcb) * erf((r + sigma) / sqrtDt4));

    const Real result(
        (term1 * (term2 + term3 + term4 + term5) + term6 + term7 + term8) / 6);
    return result;
}

struct Igbd_r_3d_params
{
    const Real sigma;
    const Real t;
    const Real D;
    const Real target;
};

static Real Igbd_r_3d_F(Real r, const Igbd_r_3d_params* params)
{
    return Igbd_r_3d(r, params->sigma, params->t, params->D) - params->target;
}

Real random_ipv_length_3d(
    RandomNumberGenerator& rng, Real const& sigma, Real const& t, Real const& D)
{
    const Real epsabs(1e-18), epsrel(1e-12);

    const Real ptot(Igbd_3d(sigma, t, D));

    Igbd_r_3d_params params = {sigma, t, D, rng.uniform(0, 1) * ptot};
    gsl_function F = {
        reinterpret_cast<typeof(F.function)>(&Igbd_r_3d_F), &params};

    Real low(sigma), high(sigma + 10 * std::sqrt(6 * D * t));

    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(gsl_root_fsolver_brent));
    gsl_root_fsolver_set(solver, &F, low, high);

    const unsigned int max_num_iter(100);
    unsigned int i(0);
    while (true)
    {
        gsl_root_fsolver_iterate(solver);

        low = gsl_root_fsolver_x_lower(solver);
        high = gsl_root_fsolver_x_upper(solver);
        int status(gsl_root_test_interval(low, high, epsabs, epsrel));

        if (status == GSL_CONTINUE)
        {
            if (i >= max_num_iter)
            {
                gsl_root_fsolver_free(solver);
                throw std::runtime_error("failed to converge");
            }
        }
        else
        {
            break;
        }

        ++i;
    }

    gsl_root_fsolver_free(solver);
    return low;
}

Position3 random_ipv_3d(
    RandomNumberGenerator& rng, Real const& sigma, Real const& t, Real const& D)
{
    const Real r(random_ipv_length_3d(rng, sigma, t, D));
    return random_spherical_uniform(rng, r);
}

} // bd

} // ecell4

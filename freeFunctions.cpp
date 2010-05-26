#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "compat.h"

#include <algorithm>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include "Logger.hpp"
#include "freeFunctions.hpp"

static Logger& _log(Logger::get_logger("freeFunctions"));

/**
   Calculates std::exp(x^2) * erfc(x)

   See asymptotic expansion here:
   http://en.wikipedia.org/wiki/Error_function
*/  
Real expxsq_erfc(Real x)
{
    Real result;

    const Real xsq(x * x);
    if(x > 26.0)
    {
        const Real M_1_SQRTPI(M_2_SQRTPI * 0.5); 

        const Real x2sq_r(1.0 / (2.0 * xsq));  // 2 / (2 x)^2

        /*
          up to second term in the expansion.
          abs err ~= 9e-8 at x == 20, 3e-8 at x == 25

          the third term 
          - (8 / (x2sq * x2sq * x2sq))       
          and beyond doesn't have a major contribution for large x.
        */

        result = (M_1_SQRTPI / x) * 
            (1.0 - x2sq_r +      // term 1
              x2sq_r * x2sq_r);  // term 2
    }
    else
    {
        result = std::exp(xsq) * erfc(x);
    }

    return result;
}


/**
   W(a, b) := std::exp(2 a b + b^2) erfc(a + b)
*/
Real W(Real a, Real b)
{
    // std::exp(2 a b + b^2) erfc(a + b) == 
    //               std::exp(- a^2) std::exp((a + b)^2) erfc(a + b)
    return std::exp(- a * a) * expxsq_erfc(a + b);
}

Real __p_irr(Real r, Real t, Real r0, Real kf, Real D, Real sigma, Real alpha)
{
    //  printf("irrp %g %g %g\n",r,r0,t);
    const Real sqrtD(std::sqrt(D));

    const Real Dt4(4.0 * D * t);
    const Real r_plus_r0_minus_2sigma(r + r0 - 2.0 * sigma);

    const Real num1(std::exp(- gsl_pow_2(r - r0) / Dt4));
    const Real num2(std::exp(- gsl_pow_2(r_plus_r0_minus_2sigma) / Dt4));
    const Real num3(W(r_plus_r0_minus_2sigma / std::sqrt(Dt4), 
                        alpha * std::sqrt(t)));

    const Real num((num1 + num2) / std::sqrt(4.0 * M_PI * t) -  alpha * num3);

    const Real den(4.0 * M_PI * r * r0 * sqrtD);

    const Real result(num / den);

    const Real jacobian(4.0 * M_PI * r * r);

    return result * jacobian;
}

Real p_irr(Real r, Real t, Real r0, Real kf, Real D, Real sigma)
{
    const Real kD(4.0 * M_PI * sigma * D);
    const Real alpha((1.0 + (kf / kD)) * (std::sqrt(D) / sigma));

    const Real p(__p_irr(r, t, r0, kf, D, sigma, alpha));

    return p;
}


Real p_survival_irr(Real t, Real r0, Real kf, Real D, Real sigma)
{
    const Real kD(4.0 * M_PI * sigma * D);
    const Real alpha((1.0 + (kf / kD)) * (std::sqrt(D) / sigma));

    const Real p(__p_reaction_irr(t, r0, kf, D, sigma, alpha, kD));

    return 1.0 - p;
}

Real __p_reaction_irr(Real t, Real r0, Real kf, Real D, Real sigma,
                       Real alpha, Real kD)
{
    const Real sqrtt(std::sqrt(t));
    const Real sqrtD(std::sqrt(D));

    const Real r0_m_sigma_over_sqrt4D_t((r0 - sigma) 
                                         / ((sqrtD + sqrtD) * sqrtt));

    const Real Wf(W(r0_m_sigma_over_sqrt4D_t, alpha * sqrtt));
    const Real factor(sigma * kf / (r0 * (kf + kD)));

    return factor * (erfc(r0_m_sigma_over_sqrt4D_t) - Wf);
}


Real 
__p_reaction_irr_t_inf(Real r0, Real kf, Real sigma, Real kD)
{
    const Real kf_kD_r0((kf + kD) * r0);
    return 1 - (kf_kD_r0 - kf * sigma) / kf_kD_r0;
}


Real p_survival_nocollision(Real t, Real r0, Real D, Real a)
{
    const Real Dt(D * t);
    const Real asq(a * a);
    const Real a_r(1.0 / a);
    const Real asq_r(a_r * a_r);

    const Real PIr0(M_PI * r0);

    const Real angle_factor(PIr0 * a_r);
    const Real exp_factor(- Dt * M_PI * M_PI * asq_r);

    const Real TOLERANCE(1e-8);

    const unsigned int i_max(
        std::max(static_cast<unsigned int>(
                      std::ceil(std::sqrt(Dt * M_PI * M_PI 
                                  + asq * std::log(1.0 / TOLERANCE) / Dt) *
                            M_1_PI)), 2u));

    Real p(0.0);
    Real sign(1.0);
    unsigned int i(1);
    while(true)
    {
        const Real term(sign * 
                         std::exp(exp_factor * i * i) * 
                         std::sin(angle_factor * i) / i);
        
        p += term;

        if(i >= i_max)
        {
            break;
        }

        sign = -sign;
        ++i;
    }

    const Real factor((a + a) / PIr0);

    return p * factor;
}

Real dp_survival_nocollision(Real t, Real r0, Real D, Real a)
{
    const Real Dt(D * t);
    const Real asq(a * a);
    const Real a_r(1.0 / a);
    const Real asq_r(a_r * a_r);

    const Real PIr0(M_PI * r0);

    const Real angle_factor(PIr0 * a_r);
    const Real exp_factor(- Dt * M_PI * M_PI * asq_r);

    const Real TOLERANCE(1e-8);

    const unsigned int i_max(
        std::max(static_cast<unsigned int>(
                      std::ceil(std::sqrt(Dt * M_PI * M_PI 
                                  + asq * std::log(1.0 / TOLERANCE) / Dt) *
                            M_1_PI)), 2u));

    Real p(0.0);
    Real sign(- 1.0);
    unsigned int i(1);
    while(true)
    {
        const Real term(sign * 
                         std::exp(exp_factor * i * i) * 
                         std::sin(angle_factor * i) * i);
        
        p += term;

        if(i >= i_max)
        {
            break;
        }

        sign = -sign;
        ++i;
    }

    const Real factor(D * (M_PI + M_PI) / (a * r0));

    return p * factor;
}

Real p_theta_free(Real theta, Real r, Real r0, Real t, Real D)
{
    Real sin_theta;
    Real cos_theta;
    sincos(theta, &sin_theta, &cos_theta);

    const Real Dt4(4.0 * D * t);
    const Real Dt4Pi(Dt4 * M_PI);

    const Real term1(std::exp(- (r * r - 2.0 * cos_theta * r * r0 + r0 * r0) / 
                           Dt4));
    const Real term2(1.0 / std::sqrt(Dt4Pi * Dt4Pi * Dt4Pi));

    return term1 * term2 * sin_theta; // jacobian
}

Real ip_theta_free(Real theta, Real r, Real r0, Real t, Real D)
{
    const Real Dt(D * t);
    const Real Dt2(Dt + Dt);
    const Real rr0(r * r0);

    const Real rr0_over_2Dt(rr0 / Dt2);

    const Real rsqr0sq_over_4Dt((r * r + r0 * r0) / (Dt2 + Dt2));

    const Real term1(expm1(rr0_over_2Dt 
                             - rsqr0sq_over_4Dt));
    const Real term2(expm1(rr0_over_2Dt * cos(theta) 
                             - rsqr0sq_over_4Dt));

    const Real den(4.0 * std::sqrt(M_PI * M_PI * M_PI * Dt) * rr0);

    return (term1 - term2) / den;
}

Real g_bd(Real r, Real sigma, Real t, Real D)
{
    const Real Dt4(4.0 * D * t);
    const Real mDt4_r(- 1.0 / Dt4);
    const Real sqrtDt4(std::sqrt(Dt4));
    const Real sqrtDt4_r(1.0 / sqrtDt4);
    const Real sqrtPi(std::sqrt(M_PI));

    const Real rps(r + sigma);
    const Real rms(r - sigma);

    const Real term1((std::exp(rps * rps * mDt4_r) - 
                        std::exp(rms * rms * mDt4_r)) * sqrtDt4 / 
                      (sqrtPi * r));
    const Real term2(erf(rps * sqrtDt4_r) - erf(rms * sqrtDt4_r));

    return 0.5 * (term1 + term2) * r * r;
}
    
Real I_bd(Real sigma, Real t, Real D)
{
    const Real sqrtPi(std::sqrt(M_PI));

    const Real Dt(D * t);
    const Real Dt2(Dt + Dt);
    const Real sqrtDt(std::sqrt(Dt));
    const Real sigmasq(sigma * sigma);

    const Real term1(1.0 / (3.0 * sqrtPi));
    const Real term2(sigmasq - Dt2);
    const Real term3(Dt2 - 3.0 * sigmasq);
    const Real term4(sqrtPi * sigmasq * sigma * erfc(sigma / sqrtDt));

    const Real result(term1 * (- sqrtDt *
                                 (term2 * std::exp(- sigmasq / Dt) + term3)
                                 + term4));
    
    return result;
}


Real I_bd_r(Real r, Real sigma, Real t, Real D)
{
    const Real sqrtPi(std::sqrt(M_PI));

    const Real Dt(D * t);
    const Real Dt2(Dt + Dt);
    const Real Dt4(Dt2 + Dt2);
    const Real sqrtDt(std::sqrt(Dt));
    const Real sqrtDt4(std::sqrt(Dt4));
    const Real sigmasq(sigma * sigma);

    const Real sigmacb(sigmasq * sigma);
    const Real rcb(gsl_pow_3(r));

    const Real rsigma(r * sigma);

    const Real rps_sq(gsl_pow_2(r + sigma));
    const Real rms_sq(gsl_pow_2(r - sigma));

    const Real term1(- 2.0 * sqrtDt / sqrtPi);
    const Real term2(std::exp(- sigmasq / Dt) * (sigmasq - Dt2));
    const Real term3(- std::exp(- rps_sq / Dt4) * (rms_sq + rsigma - Dt2));
    const Real term4(std::exp(- rms_sq / Dt4) * (rps_sq - rsigma - Dt2));
    const Real term5(- sigmasq * 3.0 + Dt2);

    const Real term6((sigmacb - rcb) * erf((r - sigma) / sqrtDt4));
    const Real term7(- (sigmacb + sigmacb) * erf(sigma / sqrtDt));
    const Real term8((sigmacb + rcb) * erf((r + sigma) / sqrtDt4));

    const Real result((term1 * (term2 + term3 + term4 + term5)
                         // + sigmasq + rsigma + rsigma - Dt2)//expm1
                         + term6 + term7 + term8) / 6.0);
    
    return result;
}


struct g_bd_params
{ 
    const Real sigma;
    const Real t;
    const Real D;
    const Real target;
};


static Real I_gbd_r_F(Real r, const g_bd_params* params)
{
    const Real sigma(params->sigma);
    const Real t(params->t);
    const Real D(params->sigma);
    const Real target(params->target);

    return I_bd_r(r, sigma, t, D) - target;
}

Real drawR_gbd(Real rnd, Real sigma, Real t, Real D)
{
    const Real I(I_bd(sigma, t, D));

    g_bd_params params = { sigma, t, D, rnd * I };

    gsl_function F =
    {
        reinterpret_cast<typeof(F.function)>(&I_gbd_r_F),
        &params
    };

    Real low(sigma);
    Real high(sigma + 10.0 * std::sqrt (6.0 * D * t));

    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));
    gsl_root_fsolver_set(solver, &F, low, high);

    const unsigned int maxIter(100);

    unsigned int i(0);
    while(true)
    {
        gsl_root_fsolver_iterate(solver);

        low = gsl_root_fsolver_x_lower(solver);
        high = gsl_root_fsolver_x_upper(solver);
        int status(gsl_root_test_interval(low, high, 1e-18, 1e-12));

        if(status == GSL_CONTINUE)
        {
            if(i >= maxIter)
            {
                gsl_root_fsolver_free(solver);
                _log.error("drawR_gbd: failed to converge");
                throw std::exception();
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

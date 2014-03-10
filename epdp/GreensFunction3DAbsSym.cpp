#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "compat.h"

#include <sstream>
#include <exception>
#include <vector>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_roots.h>

#include "findRoot.hpp"
#include "GreensFunction3DAbsSym.hpp"

/**
  EllipticTheta[4,0,q]

  Efficiently calculate EllipticTheta[4,0,q] for q < 1.0.
*/

Real GreensFunction3DAbsSym::ellipticTheta4Zero(Real q)
{
    if (fabs(q) > 1.0)
    {
        throw std::invalid_argument((boost::format("fabs(%.16g) <= 1.0") % q).str());
    }

    // et4z(1 - 1e4) ~= 7.2e-23
    // et4z(1e-15) ~= 1 - 2e-15
    // et4z(1e-16) ~= 1 - 2.2e-16
    // et4z(1e-17) ~= 1 - (zero)

    const Integer N(1000);
    Real value(1.0);
    
    Real q_n(q);
    Real q_2n(1.0);
    
    for (Integer n(1); n <= N; ++n)
    {
        const Real term2(1.0 - q_2n * q);  // q^(2n-1) = (q^(n-1))^2 * q
        
        q_2n = q_n * q_n;
        
        const Real term1(1.0 - q_2n); // q^2n
        
        const Real term(term1 * term2 * term2);
        const Real value_prev(value);
        value *= term;
        
        // here only absolute error is checked because it is good enough
        // for our use.  (it's compared with 1 in p_survival).
        if (fabs(value - value_prev) < 1e-18) 
        {
            // normal exit.
            return value;
        }
        
        q_n *= q;  // q_(++n)
    }
    
    log_.warn("ellipticTheta4Zero: didn't converge");
    return value;
}


Real GreensFunction3DAbsSym::p_survival(Real t) const
{
    const Real D(getD());
    const Real a(geta());
    const Real asq(a * a);
    const Real PIsq(M_PI * M_PI);

    const Real q(- D * PIsq * t / asq);
    return 1.0 - ellipticTheta4Zero(exp(q));
} 

Real GreensFunction3DAbsSym::p_int_r_free(Real r, Real t) const
{
    const Real D(getD());
    const Real Dt(D * t);
    const Real sqrtDt(sqrt(Dt));
    const Real sqrtPI(sqrt(M_PI));

    return erf(r / (sqrtDt + sqrtDt))
        - r * exp(- r * r / (4.0 * Dt)) / (sqrtPI * sqrtDt);
}

Real GreensFunction3DAbsSym::p_int_r(Real r, Real t) const
{
    Real value(0.0);

    const Real a(geta());
    const Real p_free(this->p_int_r_free(r, t));

    // p_int_r is always smaller than p_free.
    if (fabs(p_free) < CUTOFF)
    {
        return 0.0;
    }

    const Real D(getD());
    const Real asq(a * a);
    const Real PIsq(M_PI * M_PI);

    const Real PIr(M_PI * r);
    const Real PIr_a(PIr / a);
    const Real DtPIsq_asq(D * t * PIsq / asq);
    
    const Real factor(2.0 / (a * M_PI));

    const Real maxn((a / M_PI) * sqrt(log(exp(DtPIsq_asq) / CUTOFF) / 
                                          (D * t)));

    const Integer N_MAX(10000);

    const Integer N(std::min(static_cast<Integer>(ceil(maxn) + 1),
                               N_MAX));
    if (N == N_MAX)
    {
        log_.warn("p_int_r: didn't converge");
    }
    

    for (Integer n(1); n <= N; ++n)
    {
        const Real term1(exp(- n * n * DtPIsq_asq));
      
        const Real angle_n(n * PIr_a);
        Real sin_n;
        Real cos_n;
        sincos(angle_n, &sin_n, &cos_n);
        const Real term2(a * sin_n);
        const Real term3(n * PIr * cos_n);

        const Real term(term1 * (term2 - term3) / n);
        value += term;
    }

    return value * factor;
} 

Real GreensFunction3DAbsSym::p_r_fourier(Real r, Real t) const 
{
    Real value(0.0);

    const Real D(getD());
    const Real a(geta());
    const Real asq(a * a);
    const Real PIsq(M_PI * M_PI);

    const Integer N(100);

    long int n(1);
    for (;;)
    {
        const Real term1(exp(- (PIsq * r * r + asq * n*n) / 
                               (4.0 * D * PIsq * t)));

        const Real term2(M_PI * r * 
                          exp(gsl_sf_lncosh(a * r * n / 
                                              (2.0 * D * M_PI * t))));

        const Real term3(a * n *
                          exp(gsl_sf_lnsinh(a * r * n / 
                                              (2.0 * D * M_PI * t))));


        const Real term(term1 * r * (term2 - term3));
        value += term;

        if (fabs(value) * 1e-8 > fabs(term))
        {
            break;
        }

        if (n > N)
        {
            log_.warn("p_r_fourier: didn't converge; n = %d, value = %.16g", n, value);
            break;
        }

        ++n;
    }

    const Real factor(1.0 / (sqrt(2) * PIsq * pow(D * t, 1.5)));

    return value * factor;
} 

struct p_survival_params
{
    const GreensFunction3DAbsSym* const gf;
    const Real rnd;
};

static Real p_survival_F(Real t, p_survival_params const* params)
{
    return params->rnd - params->gf->p_survival(t);
}


Real GreensFunction3DAbsSym::drawTime(Real rnd) const
{
    const Real D(getD());

    if (rnd >= 1.0 || rnd < 0.0)
    {
        throw std::invalid_argument((boost::format("0.0 <= %.16g < 1.0") % rnd).str());
    }

    const Real a(geta());

    if (D == 0.0 || a == INFINITY)
    {
        return INFINITY;
    }

    if (a == 0.0)
    {
        return 0.0;
    }

    p_survival_params params = { this, rnd };

    gsl_function F = 
        {
            reinterpret_cast<typeof(F.function)>(&p_survival_F),
            &params 
        };

    const Real t_guess(a * a / (6. * D));

    Real low(t_guess);
    Real high(t_guess);

    const Real value(GSL_FN_EVAL(&F, t_guess));

    if (value < 0.0)
    {
        high *= 10;

        for (;;)
        {
            const Real high_value(GSL_FN_EVAL(&F, high));
            
            if (high_value >= 0.0)
            {
                break;
            }

            if (fabs(high) >= t_guess * 1e6)
            {
                throw std::runtime_error(
                    (boost::format("couldn't adjust high. F(%.16g) = %.16g; %s") %
                       high % GSL_FN_EVAL(&F, high) %
                       boost::lexical_cast<std::string>(*this)).str());
            }
            high *= 10;
        }
    }
    else
    {
        Real low_value_prev(value);
        low *= .1;

        for (;;)
        {
            const Real low_value(GSL_FN_EVAL(&F, low));
            
            if (low_value <= 0.0)
            {
                break;
            }
            
            if (fabs(low) <= t_guess * 1e-6 ||
                fabs(low_value - low_value_prev) < CUTOFF)
            {
                log_.info("couldn't adjust high. F(%.16g) = %.16g; %s",
                          low, GSL_FN_EVAL(&F, low),
                          boost::lexical_cast<std::string>(*this).c_str());
                log_.info("returning low (%.16g)", low);
                return low;
            }
            low_value_prev = low_value;
            low *= .1;
        }
    }


    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));

    const Real t(findRoot(F, solver, low, high, 1e-18, 1e-12,
                            "GreensFunction3DAbsSym::drawTime"));

    gsl_root_fsolver_free(solver);

    return t;
}

struct p_r_params
{
    const GreensFunction3DAbsSym* const gf;
    const Real t;
    const Real target;
};

static Real p_r_free_F(Real r, p_r_params const* params)
{
    return params->gf->p_int_r_free(r, params->t) - params->target;
}


static Real p_r_F(Real r, p_r_params const* params)
{
    return params->gf->p_int_r(r, params->t) - params->target;
}

Real GreensFunction3DAbsSym::drawR(Real rnd, Real t) const 
{
    if (rnd >= 1.0 || rnd < 0.0)
    {
        throw std::invalid_argument((boost::format("0.0 <= %.16g < 1.0") % rnd).str());
    }

    if (t < 0.0)
    {
        throw std::invalid_argument((boost::format("%.16g < 0.0") % t).str());
    }

    const Real a(geta());
    const Real D(getD());

    if (a == 0.0 || t == 0.0 || D == 0.0)
    {
        return 0.0;
    }

    const Real thresholdDistance(this->CUTOFF_H * sqrt(6.0 * D * t));

    gsl_function F;
    Real psurv;

    if (a <= thresholdDistance)
    {
        //psurv = p_survival(t);  // this causes a problem when p_survival is very small.
        psurv = p_int_r(a, t);

        if (psurv == 0.0)
        {
            return a;
        }

        assert(psurv >= 0.0);

        F.function = reinterpret_cast<typeof(F.function)>(&p_r_F);
    }
    else
    {
        // p_int_r < p_int_r_free
        if (p_int_r_free(a, t) < rnd)
        {
            log_.info("p_int_r_free(a, t) < rnd, returning a");
            return a;
        }

        psurv = 1.0;
        F.function = reinterpret_cast<typeof(F.function)>(&p_r_free_F);
    }

    const Real target(psurv * rnd);
    p_r_params params = { this, t, target };

    F.params = &params;

    const Real low(0.0);
    const Real high(a);
    //const Real high(std::min(thresholdDistance, a));

    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));

    const Real r(findRoot(F, solver, low, high, 1e-18, 1e-12,
                            "GreensFunction3DAbsSym::drawR"));
  
    gsl_root_fsolver_free(solver);

    return r;
}

std::string GreensFunction3DAbsSym::dump() const
{
    return (boost::format("D=%.16g, a=%.16g") % getD() % geta()).str();
}

Logger& GreensFunction3DAbsSym::log_(
        Logger::get_logger("GreensFunction3DAbsSym"));

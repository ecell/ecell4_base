#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "compat.h"

#include <stdexcept>
#include <vector>
#include <sstream>
#include <algorithm>

#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_lambert.h>

#include "funcSum.hpp"
#include "freeFunctions.hpp"
#include "SphericalBesselGenerator.hpp"
#include "GreensFunction3DAbs.hpp"

typedef GreensFunction3DAbs GF3DA;

GF3DA::GreensFunction3DAbs(Real D, Real r0, Real a) 
    : GreensFunction3DRadAbsBase(D, 0., r0, 0.), a(a)
{
    if (a < 0.0)
    {
        throw std::invalid_argument((boost::format("a >= 0.0 : a=%.16g") % a).str());
    }
}

GF3DA::~GreensFunction3DAbs()
{
    ; // do nothing
}

Real
GF3DA::p_survival(Real t) const
{
    const Real D(getD());
    const Real a(geta());

    return p_survival_nocollision(t, r0, D, a);
}


Real
GF3DA::dp_survival(Real t) const
{
    const Real D(getD());
    const Real a(geta());

    return dp_survival_nocollision(t, r0, D, a);
}

Real
GF3DA::p_int_r(Real r, Real t) const
{
    const Real D(getD());
    const Real a(geta());

    const Real Dt(D * t);
    const Real asq(a * a);
    const Real a_r(1.0 / a);
    const Real asq_r(a_r * a_r);

    const Real PIr0(M_PI * r0);
    const Real PIr(M_PI * r);

    const Real r0_angle_factor(PIr0 * a_r);
    const Real r_angle_factor(PIr * a_r);
    const Real exp_factor(- Dt * M_PI * M_PI * asq_r);

    const unsigned int i_max(
        std::max(static_cast<unsigned int>(
                      ceil(sqrt( 1.0 - asq / M_PI / M_PI
                                    * log(TOLERANCE) / Dt ) 
                          )), 
                  2u));

    Real p(0.0);
    unsigned int i(1);
    for (;;)
    {
        Real sin_r;
        Real cos_r;
        sincos(r_angle_factor * i, &sin_r, &cos_r);

        const Real isq(i * i);

        const Real term1(exp(exp_factor * isq) * 
                          sin(r0_angle_factor * i));
        const Real term2(a * sin_r - PIr * i * cos_r);
        const Real term(term1 * term2 / isq);
        
        p += term;

        if (i >= i_max)
        {
            break;
        }

        ++i;
    }

    const Real factor(M_2_PI / PIr0);

    return p * factor;
}

struct p_survival_params
{ 
    const GF3DA* const gf;
    const Real rnd;
};

Real
static p_survival_F(Real t, p_survival_params const* params)
{
    const GF3DA* const gf(params->gf); 
    const Real rnd(params->rnd);

    return rnd - gf->p_survival(t);
}

struct p_int_r_params
{ 
    const GF3DA* const gf;
    const Real t;
    const Real rnd;
};

static Real p_int_r_F(Real r, p_int_r_params const* params)
{
    const GF3DA* const gf(params->gf); 
    const Real t(params->t);
    const Real rnd(params->rnd);

    return gf->p_int_r(r, t) - rnd;
}

Real 
GF3DA::p_n_alpha(unsigned int i, unsigned int n, Real r, Real t) const
{
    const Real a(geta());

    const Real mDt(- getD() * t);

    // j = a alpha -> alpha = j / a
    const Real aalpha(gsl_sf_bessel_zero_Jnu(static_cast<Real>(n) + 0.5, 
                                               i + 1));
    const Real alpha(aalpha / a);

    const Real term1(exp(mDt * alpha * alpha));

    const SphericalBesselGenerator& s(SphericalBesselGenerator::instance());

    const Real jr( s.j(n,   r * alpha));
    const Real jr0(s.j(n,   r0 * alpha));
    const Real ja2(s.j(n+1,   aalpha));

    const Real num(jr * jr0);
    const Real den(ja2 * ja2);

    const Real result(term1 * num / den);

    return result;
}


Real 
GF3DA::p_n(Integer n, Real r, Real t) const
{
    const Real p(funcSum(boost::bind(&GF3DA::p_n_alpha, this, _1, n, r, t),
                         MAX_ALPHA_SEQ));

    return p;
}

void
GF3DA::makep_nTable(RealVector& p_nTable, Real r, Real t) const
{
    const Real a(geta());

    p_nTable.clear();

    const Real factor(1.0 / (2.0 * M_PI * gsl_pow_3(a))); 

    const Real p_0(p_n(0, r, t) * factor);
    p_nTable.push_back(p_0);

    if (p_0 == 0)
    {
        return;
    }

    const Real threshold(fabs(p_0 * THETA_TOLERANCE * 1e-1 ));

    Real p_n_prev_abs(fabs(p_0));
    unsigned int n(1);
    for (;;) 
    {
        Real p_n(this->p_n(n, r, t) * factor);

        if (! isfinite(p_n))
        {
            log_.error("makep_nTable: invalid value: %.16g (n=%d)", p_n, n);
            break;
        }

        p_nTable.push_back(p_n);

        const Real p_n_abs(fabs(p_n));
        // truncate when converged enough.
        if (p_n_abs <= threshold &&
            p_n_prev_abs <= threshold &&
            p_n_abs <= p_n_prev_abs)
        {
            break;
        }
        

        if (n >= MAX_ORDER)
        {
            break;
        }
        
        ++n;
        p_n_prev_abs = p_n_abs;
    }

}



static Real
p_theta_i(unsigned int n,
          GF3DA::RealVector const& p_nTable,
          GF3DA::RealVector const& lgndTable)
{
    return p_nTable[n] * lgndTable[n] * (2 * n + 1);
}

Real GF3DA::p_theta_table(Real theta, Real r, Real t, RealVector const& p_nTable) const
{
    const unsigned int tableSize(p_nTable.size());

    RealVector lgndTable(tableSize);

    Real sin_theta;
    Real cos_theta;
    sincos(theta, &sin_theta, &cos_theta);
    gsl_sf_legendre_Pl_array(tableSize-1, cos_theta, &lgndTable[0]);

    return funcSum_all(
        boost::bind(&p_theta_i, _1, p_nTable, lgndTable),
        tableSize) * sin_theta;
}


Real GF3DA::p_theta(Real theta, Real r, Real t) const 
{
    {
        const Real a(geta());
        
        if (!(theta >= 0.0 && theta <= M_PI))
        {
            throw std::invalid_argument((boost::format("theta >= 0.0 && theta <= M_PI : theta=%.16g, M_PI=%.16g") % theta % M_PI).str());
        }

        if (!(r >= 0 && r < a))
        {
            throw std::invalid_argument((boost::format("r >= 0 && r < a : r=%.16g, a=%.16g") % r % a).str());
        }

        if (!(r0 >= 0 && r0 < a))
        {
            throw std::invalid_argument((boost::format("r0 >= 0 && r0 < a : r0=%.16g, a=%.16g") % r0 % a).str());
        }

        if (!(t >= 0.0))
        {
            throw std::invalid_argument((boost::format("t >= 0.0 : t=%.16g") % t).str());
        }

    }

    if (t == 0.0)
    {
        return 0.0;
    }

    
    RealVector p_nTable;

    makep_nTable(p_nTable, r, t);

    const Real p(p_theta_table(theta, r, t, p_nTable));

    return p;
}



Real GF3DA::ip_theta(Real theta, Real r, Real t) const
{
    {
        const Real a(geta());
        
        if (!(theta >= 0.0 && theta <= M_PI))
        {
            throw std::invalid_argument((boost::format("theta >= 0.0 && theta <= M_PI : theta=%.16g, M_PI=%.16g") % theta % M_PI).str());
        }

        // r \in (sigma, a)
        if (!(r >= 0.0 && r < a))
        {
            throw std::invalid_argument((boost::format("r >= 0.0 && r < a : r=%.16g, a=%.16g") % r % a).str());
        }

        if (!(r0 >= 0.0 && r0 < a))
        {
            throw std::invalid_argument((boost::format("r0 >= 0.0 && r0 < a : r0=%.16g, a=%.16g") % r0 % a).str());
        }

        if (!(t >= 0.0))
        {
            throw std::invalid_argument((boost::format("t >= 0.0 : t=%.16g") % t).str());
        }

    }

    if (t == 0.0 || theta == 0.0)
    {
        return 0.0;
    }

    RealVector p_nTable;

    makep_nTable(p_nTable, r, t);

    const Real p(ip_theta_table(theta, r, t, p_nTable));

    return p;
}

static Real
ip_theta_i(unsigned int n,
           GF3DA::RealVector const& p_nTable,
           GF3DA::RealVector const& lgndTable1)
{
    // lgndTable1 is offset by 1; lgndTable1[0] is for n=-1.

    const Real lgnd_n_m1(lgndTable1[n]);   // n-1
    const Real lgnd_n_p1(lgndTable1[n+2]); // n+1
    
    return p_nTable[n] * (lgnd_n_m1 - lgnd_n_p1);// / (1.0 + 2.0 * n);
}


Real 
GF3DA::ip_theta_table(
    Real theta, Real r, Real t, RealVector const& p_nTable) const
{
    const unsigned int tableSize(p_nTable.size());

    RealVector pTable;
    pTable.reserve(tableSize);

    const Real cos_theta(cos(theta));

    // LgndTable is offset by 1 to incorporate the n=-1 case.
    // For ex: LgndTable[0] is for n=-1, lgndTable[1] is for n=0 ...

    RealVector lgndTable1(tableSize + 2);
    lgndTable1[0] = 1.0;  // n = -1
    gsl_sf_legendre_Pl_array(tableSize, cos_theta, &lgndTable1[1]);


    return funcSum_all(
        boost::bind(&ip_theta_i, _1, p_nTable, lgndTable1),
        tableSize);
}

struct GF3DA::ip_theta_params
{ 
    GF3DA const* const gf;
    const Real r;
    const Real t;
    RealVector const& p_nTable;
    const Real value;
};

Real GF3DA::ip_theta_F(Real theta, ip_theta_params const* params)
{
    const GF3DA* const gf(params->gf); 
    const Real r(params->r);
    const Real t(params->t);
    const RealVector& p_nTable(params->p_nTable);
    const Real value(params->value);

    return gf->ip_theta_table(theta, r, t, p_nTable) - value;
}


Real 
GF3DA::dp_n_alpha(unsigned int i, unsigned int n, Real t) const
{
    const Real a(geta());

    const Real mDt(- getD() * t);

    const Real 
        aalpha(gsl_sf_bessel_zero_Jnu(static_cast<Real>(n) + 0.5, i + 1));
    const Real alpha(aalpha / a);

    const Real term1(exp(mDt * alpha * alpha) * alpha);

    const SphericalBesselGenerator& s(SphericalBesselGenerator::instance());

    const Real jr0(s.j(n,   r0 * alpha));
    const Real ja2(s.j(n+1,   aalpha));

    const Real result(term1 * jr0 / ja2);

    return result;
}



Real 
GF3DA::dp_n(Integer n, Real t) const
{
    const Real 
        p(funcSum(boost::bind(&GF3DA::dp_n_alpha, this, _1, n, t),
                  MAX_ALPHA_SEQ));

    return p;
}


void GF3DA::makedp_nTable(RealVector& p_nTable, Real t) const
{
    p_nTable.clear();

    const Real factor(- getD() / (2.0 * M_PI * gsl_pow_3(a)));

    const Real p_0(dp_n(0, t) * factor);
    p_nTable.push_back(p_0);

    if (p_0 == 0)
    {
        return;
    }

    const Real threshold(fabs(THETA_TOLERANCE * p_0 * 1e-1));

    Real p_n_prev_abs(fabs(p_0));
    unsigned int n(1);
    for (;;)
    {
        Real p_n(dp_n(n, t) * factor);

        if (! isfinite(p_n))
        {
            log_.error("makedp_nTable: invalid value: %.16g (n=%d)", p_n, n);
            break;
        }

        p_nTable.push_back(p_n);

        const Real p_n_abs(fabs(p_n));
        // truncate when converged enough.
        if (p_n_abs <= threshold &&
            p_n_prev_abs <= threshold &&
            p_n_abs <= p_n_prev_abs)
        {
            break;
        }
        

        if (n >= MAX_ORDER)
        {
            break;
        }
        
        ++n;
        p_n_prev_abs = p_n_abs;
    }

}

Real 
GF3DA::dp_theta(Real theta, Real r, Real t) const 
{
    {
        const Real a(geta());
        
        if (!(theta >= 0.0 && theta <= M_PI))
        {
            throw std::invalid_argument((boost::format("theta >= 0.0 && theta <= M_PI : theta=%.16g, M_PI=%.16g") % theta % M_PI).str());
        }


        // r \in [ sigma, a ]  ;  unlike p_theta,
        // defined at r == sigma and r == a.
        if (!(r >= 0.0 && r <= a))
        {
            throw std::invalid_argument((boost::format("r >= 0.0 && r <= a : r=%.16g, a=%.16g") % r % a).str());
        }

        if (!(r0 >= 0.0 && r0 < a))
        {
            throw std::invalid_argument((boost::format("r0 >= 0.0 && r0 < a : r0=%.16g, a=%.16g") % r0 % a).str());
        }

        if (!(t >= 0.0))
        {
            throw std::invalid_argument((boost::format("t >= 0.0 : t=%.16g") % t).str());
        }

    }

    if (t == 0.0)
    {
        return 0.0;
    }

    RealVector p_nTable;

    makedp_nTable(p_nTable, t);

    const Real p(p_theta_table(theta, r, t, p_nTable));

    return p;
}

Real 
GF3DA::idp_theta(Real theta, Real r, Real t) const
{
    {
        const Real a(geta());
        
        if (!(theta >= 0.0 && theta <= M_PI))
        {
            throw std::invalid_argument((boost::format("theta >= 0.0 && theta <= M_PI : theta=%.16g, M_PI=%.16g") % theta % M_PI).str());
        }

        // r \in [ sigma, a ]
        if (!(r >= 0.0 && r <= a))
        {
            throw std::invalid_argument((boost::format("r >= 0.0 && r <= a : r=%.16g, a=%.16g") % r % a).str());
        }

        if (!(r0 >= 0.0 && r0 < a))
        {
            throw std::invalid_argument((boost::format("r0 >= 0.0 && r0 < a : r0=%.16g, a=%.16g") % r0 % a).str());
        }

        if (!(t >= 0.0))
        {
            throw std::invalid_argument((boost::format("t >= 0.0 : t=%.16g") % t).str());
        }

    }

    if (t == 0.0 || theta == 0.0)
    {
        return 0.0;
    }

    RealVector p_nTable;

    makedp_nTable(p_nTable, t);

    const Real p(ip_theta_table(theta, r, t, p_nTable));

    return p;
}

Real 
GF3DA::drawTime(Real rnd) const
{
   const Real a(geta());

   if (!(rnd <= 1.0 && rnd >= 0.0))
   {
       throw std::invalid_argument((boost::format("rnd <= 1.0 && rnd >= 0.0 : rnd=%.16g") % rnd).str());
   }

   if (!(r0 >= 0.0 && r0 <= a))
   {
       throw std::invalid_argument((boost::format("r0 >= 0.0 && r0 <= a : r0=%.16g, a=%.16g") % r0 % a).str());
   }


   if (r0 == a || a == 0.0)
   {
       return 0.0;
   }

   Real low(1e-6);
   Real high(1.0);

   p_survival_params params = { this, rnd };

   gsl_function F = {
       reinterpret_cast<typeof(F.function)>(&p_survival_F),
       &params 
   };

    // adjust high and low to make sure that f(low) and f(high) straddle.
    while(GSL_FN_EVAL(&F, high) < 0.0)
    {
        high *= 10;
        log_.info("drawTime: adjusting high: %.16g", high);
        if (fabs(high) >= 1e10)
        {
            throw std::runtime_error(
                (boost::format("couldn't adjust high. F(%.16g) = %.16g; r0=%.16g, %s") %
                    high % GSL_FN_EVAL(&F, high) % r0 % dump()).str());
        }
    }

    Real low_value(GSL_FN_EVAL(&F, low));
    while(low_value > 0.0)
    {
        low *= .1;

        const Real low_value_new(GSL_FN_EVAL(&F, low));

        log_.info("drawTime: adjusting low: %.16g, F = %.16g", low, low_value_new);

        if (fabs(low) <= MIN_T || 
            fabs(low_value - low_value_new) < TOLERANCE) 
        {
            log_.info("couldn't adjust low.  Returning %.16g as MIN_T; "
                      "F(%.16g) = %.16g; r0 = %.16g, %s", MIN_T, low,
                      GSL_FN_EVAL(&F, low), r0, dump().c_str());
            return MIN_T;
        }

        low_value = low_value_new;
    }

    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));
    gsl_root_fsolver_set(solver, &F, low, high);

    const unsigned int maxIter(100);

    unsigned int i(0);
    for (;;)
    {
        gsl_root_fsolver_iterate(solver);
        low = gsl_root_fsolver_x_lower(solver);
        high = gsl_root_fsolver_x_upper(solver);

        const int status(gsl_root_test_interval(low, high, MIN_T, 
                                                  TOLERANCE));

        if (status == GSL_CONTINUE)
        {
            if (i >= maxIter)
            {
                gsl_root_fsolver_free(solver);
                throw std::runtime_error("drawTime: failed to converge");
            }
        }
        else
        {
            break;
        }

        ++i;
    }
  
    Real t(gsl_root_fsolver_root(solver));
    gsl_root_fsolver_free(solver);

    return t;
}

Real 
GF3DA::drawR(Real rnd, Real t) const
{
    const Real a(geta());

    if (!(rnd <= 1.0 && rnd >= 0.0))
    {
        throw std::invalid_argument((boost::format("rnd <= 1.0 && rnd >= 0.0 : rnd=%.16g") % rnd).str());
    }

    if (!(r0 >= 0.0 && r0 < a))
    {
        throw std::invalid_argument((boost::format("r0 >= 0.0 && r0 < a : r0=%.16g, a=%.16g") % r0 % a).str());
    }


    if (t == 0.0)
    {
        return r0;
    }

    const Real psurv(p_survival(t));

    p_int_r_params params = { this, t, rnd * psurv };

    gsl_function F = {
        reinterpret_cast<typeof(F.function)>(&p_int_r_F),
        &params 
    };

    Real low(0.0);
    Real high(a);

//    const Real lowvalue(GSL_FN_EVAL(&F, low ));
    const Real highvalue(GSL_FN_EVAL(&F, high));

    // No initial range guess, except the negative value check below,
    // as evaluation of p_int_r in this GF seems pretty robust.

    if (highvalue < 0.0)
    {
        log_.info("drawR: highvalue < 0.0 (%.16g). returning a (%.16g)", highvalue, a);
        return a;
    }


    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));
    gsl_root_fsolver_set(solver, &F, low, high);

    const unsigned int maxIter(100);

    unsigned int i(0);
    for (;;)
    {
        gsl_root_fsolver_iterate(solver);
        low = gsl_root_fsolver_x_lower(solver);
        high = gsl_root_fsolver_x_upper(solver);
        const int status(gsl_root_test_interval(low, high, 1e-15,
                                                  TOLERANCE));

        if (status == GSL_CONTINUE)
        {
            if (i >= maxIter)
            {
                gsl_root_fsolver_free(solver);
                throw std::runtime_error("drawR: failed to converge");
            }
        }
        else
        {
            break;
        }

        ++i;
    }
  
    const Real r(gsl_root_fsolver_root(solver));
    gsl_root_fsolver_free(solver);

    return r;
}
    
Real 
GF3DA::drawTheta(Real rnd, Real r, Real t) const
{
    Real theta;

    const Real a(geta());

    // input parameter range checks.
    if (!(rnd <= 1.0 && rnd >= 0.0))
    {
        throw std::invalid_argument((boost::format("rnd <= 1.0 && rnd >= 0.0 : rnd=%.16g") % rnd).str());
    }

    if (!(r0 >= 0.0 && r0 < a))
    {
        throw std::invalid_argument((boost::format("r0 >= 0.0 && r0 < a : r0=%.16g, a=%.16g") % r0 % a).str());
    }

    if (!(r >= 0.0 && r <= a))
    {
        throw std::invalid_argument((boost::format("r >= 0.0 && r <= a : r=%.16g, a=%.16g") % r % a).str());
    }

    if (!(t >= 0.0))
    {
        throw std::invalid_argument((boost::format("t >= 0.0 : t=%.16g") % t).str());
    }


    // t == 0 means no move.
    if (t == 0.0)
    {
        return 0.0;
    }

    RealVector p_nTable;

    if (r == geta() || r < 0.0)
    {
        makedp_nTable(p_nTable, t);
    }
    else
    {
        makep_nTable(p_nTable, r, t);
    }

    // root finding with the integrand form.

    const Real ip_theta_pi(ip_theta_table(M_PI, r, t, p_nTable));

    ip_theta_params params = { this, r, t, p_nTable, rnd * ip_theta_pi };

    gsl_function F = {
        reinterpret_cast<typeof(F.function)>(&ip_theta_F),
        &params 
    };

    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));
    gsl_root_fsolver_set(solver, &F, 0.0, M_PI);

    const unsigned int maxIter(100);

    unsigned int i(0);
    for (;;)
    {
        gsl_root_fsolver_iterate(solver);
        const Real low(gsl_root_fsolver_x_lower(solver));
        const Real high(gsl_root_fsolver_x_upper(solver));
        const int status(gsl_root_test_interval(low, high, 1e-11,
                                                  THETA_TOLERANCE));

        if (status == GSL_CONTINUE)
        {
            if (i >= maxIter)
            {
                gsl_root_fsolver_free(solver);
                throw std::runtime_error("drawTheta: failed to converge");
            }
        }
        else
        {
            break;
        }

        ++i;
    }
  
    theta = gsl_root_fsolver_root(solver);
    gsl_root_fsolver_free(solver);
    
    return theta;
}

GF3DA::EventKind GF3DA::drawEventType(Real rnd, Real t) const
{
    assert(0);
}

//
// debug
//

std::string GF3DA::dump() const
{
    std::ostringstream ss;
    ss << "D = " << getD() <<
        ", a = " << geta() << std::endl;
    return ss.str();
}

Logger& GreensFunction3DAbs::log_(
    Logger::get_logger("GreensFunction3DAbs"));

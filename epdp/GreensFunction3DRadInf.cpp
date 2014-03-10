#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdexcept>
#include <vector>
#include <sstream>
#include <cmath>

#include "compat.h"

#include <boost/bind.hpp>
#include <boost/format.hpp>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_roots.h>

#include "freeFunctions.hpp"

#include "funcSum.hpp"

#include "SphericalBesselGenerator.hpp"

#include "GreensFunction3DRadInf.hpp"



GreensFunction3DRadInf::GreensFunction3DRadInf(Real D, Real kf, Real r0, Real Sigma)
    : PairGreensFunction(D, kf, r0, Sigma),
      kD(4.0 * M_PI * getSigma() * getD()),
      alpha((1.0 + (getkf() / getkD())) * (sqrt(getD()) / getSigma()))
{
    ; // do nothing
}

GreensFunction3DRadInf::~GreensFunction3DRadInf()
{
    ; // do nothing
}

Real 
GreensFunction3DRadInf::p_corr_R(Real alpha, unsigned int n, Real r, Real t) const
{
    const Real D(this->getD());
    const Real sigma(this->getSigma());
    
    const Real ks(getkf() * sigma);
    const Real realn(static_cast<Real>(n));
    const Real ks_m_n(ks - realn);

    const Real alphasq(alpha * alpha);

    const Real term1(exp(- D * t * alphasq));

    const Real sAlpha(sigma * alpha);
    const Real rAlpha(r * alpha);
    const Real r0Alpha(r0 * alpha);

    const SphericalBesselGenerator& s(SphericalBesselGenerator::instance());
    const Real js( s.j(n,     sAlpha));
    const Real ys( s.y(n,     sAlpha));
    const Real js1(s.j(n + 1, sAlpha));
    const Real ys1(s.y(n + 1, sAlpha));
    const Real jr( s.j(n,     rAlpha));
    const Real yr( s.y(n,     rAlpha));
    const Real jr0(s.j(n,     r0Alpha));
    const Real yr0(s.y(n,     r0Alpha));

    const Real R1((ks_m_n * js + sAlpha * js1));
    const Real R2((ks_m_n * ys + sAlpha * ys1));

    const Real F1R1(R1 * jr * jr0 - R1 * yr * yr0);
    const Real F2(jr0 * yr + jr * yr0);

    const Real num(2.0 * sqrt(r * r0) *
                    alphasq * R1 * (F1R1 + F2 * R2));
    const Real den(M_PI * (R1 * R1 + R2 * R2));

    const Real result(term1 * num / den);

    assert(isfinite(result));

    return result;
}

struct GreensFunction3DRadInf::p_corr_R_params 
{ 
    const GreensFunction3DRadInf* const gf;
    unsigned int n;
    const Real r;
    const Real t; 
};

Real GreensFunction3DRadInf::p_corr_R_F(Real alpha, p_corr_R_params* params)
{
    const GreensFunction3DRadInf* const gf(params->gf); 

    const unsigned int n(params->n);
    const Real r(params->r);
    const Real t(params->t);

    return gf->p_corr_R(alpha, n, r, t);
}


Real GreensFunction3DRadInf::p_corr(Real theta, Real r, Real t) const
{
    RealVector RnTable;
    makeRnTable(RnTable, r, t);

    return p_corr_table(theta, r, t, RnTable);
}

Real GreensFunction3DRadInf::ip_corr(Real theta, Real r, Real t) const
{
    RealVector RnTable;
    makeRnTable(RnTable, r, t);

    return ip_corr_table(theta, r, t, RnTable);
}


Real GreensFunction3DRadInf::p_free(Real theta, Real r, Real t) const
{
    return p_theta_free(theta, r, r0, t, getD());
}

Real GreensFunction3DRadInf::p_survival(Real t) const
{
    return 1.0 - p_reaction(t);
}


Real GreensFunction3DRadInf::p_reaction(Real t) const
{
    const Real kf(getkf());
    const Real D(getD());
    const Real sigma(getSigma());
    const Real alpha(getalpha());
    const Real kD(getkD());

    return __p_reaction_irr(t, r0, kf, D, sigma, alpha, kD);
}

struct p_reaction_params 
{ 
    const GreensFunction3DRadInf* const gf;
    const Real rnd;
};

static Real p_reaction_F(Real t, p_reaction_params* params)
{
    const GreensFunction3DRadInf* const gf(params->gf); 
    const Real kf(gf->getkf());
    const Real D(gf->getD());
    const Real sigma(gf->getSigma());
    const Real alpha(gf->getalpha());
    const Real kD(gf->getkD());

    const Real r0(gf->getr0());
    const Real rnd(params->rnd);

    return __p_reaction_irr(t, r0, kf, D, sigma, alpha, kD ) - rnd;
}


Real 
GreensFunction3DRadInf::p_int_r(Real r, Real t) const
{
    const Real kf(getkf());
    const Real D(getD());
    const Real sigma(getSigma());
    const Real alpha(getalpha());
    const Real kD(getkD());

    const Real Dt(D * t);

    const Real kf_kD(kf + kD);
    const Real Dt4(4.0 * Dt);
    const Real sqrtDt4(sqrt(Dt4));
    const Real ksigma2(2.0 * kf * sigma);
    const Real alphasqrtt(alpha * sqrt(t));

    const Real r_r0__2s___sqrtDt4((r - 2.0 * sigma + r0) / sqrtDt4);
    const Real r_r0__sqrtDt4((r - r0) / sqrtDt4);
    const Real r0_s__sqrtDt4((r0 - sigma) / sqrtDt4);

    const Real term1((expm1(- gsl_pow_2(r_r0__2s___sqrtDt4 ))
                        - expm1(- gsl_pow_2(r_r0__sqrtDt4))) * 
                        sqrt(Dt / M_PI));

    const Real erf_r_r0__2s___sqrtDt4(erf(r_r0__2s___sqrtDt4));
    const Real term2(kf_kD * r0 * erf(r_r0__sqrtDt4) 
                      + kf_kD * r0 * erf_r_r0__2s___sqrtDt4
                      + ksigma2 * 
                      (erf(r0_s__sqrtDt4) - erf_r_r0__2s___sqrtDt4));

    const Real term3(kf * sigma * W(r0_s__sqrtDt4, alphasqrtt) 
                      - (kf * r + kD * (r - sigma)) *
                      W(r_r0__2s___sqrtDt4, alphasqrtt));

    const Real result((1 / r0) * (term1 + (1 / kf_kD) * 
                                      ((0.5 * term2) + term3)));

    return result;
}

struct p_int_r_params 
{ 
    const GreensFunction3DRadInf* const gf;
    const Real t;
    const Real rnd;
};

static Real p_int_r_F(Real r, p_int_r_params* params)
{
    const GreensFunction3DRadInf* const gf(params->gf); 

    const Real t(params->t);
    const Real rnd(params->rnd);

    return gf->p_int_r(r, t) - rnd;
}


Real GreensFunction3DRadInf::drawTime(Real rnd) const
{
    const Real sigma(this->getSigma());

    if (!(rnd < 1.0 && rnd >= 0.0))
    {
        throw std::invalid_argument((boost::format("rnd < 1.0 && rnd >= 0.0 : rnd=%.16g") % rnd).str());
    }

    if (!(r0 >= sigma))
    {
        throw std::invalid_argument((boost::format("r0 >= sigma : r0=%.16g, sigma=%.16g") % r0 % sigma).str());
    }


    Real low(1e-100);
    Real high(100);

    {
        const Real maxp(p_reaction(INFINITY));

        if(rnd >= maxp)
        {
            return INFINITY;
        }
    }

    p_reaction_params params = { this, rnd };

    gsl_function F = 
        {
            reinterpret_cast<typeof(F.function)>(&p_reaction_F),
            &params 
        };

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
        int status(gsl_root_test_interval(low, high, 1e-18, 1e-12));

        if(status == GSL_CONTINUE)
        {
            if(i >= maxIter)
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
  
    const Real r(gsl_root_fsolver_root(solver));
    gsl_root_fsolver_free(solver);

    return r;
} 




Real GreensFunction3DRadInf::drawR(Real rnd, Real t) const
{
    const Real sigma(this->getSigma());
    const Real D(this->getD());

    if (!(rnd < 1.0 && rnd >= 0.0))
    {
        throw std::invalid_argument((boost::format("rnd < 1.0 && rnd >= 0.0 : rnd=%.16g") % rnd).str());
    }

    if (!(r0 >= sigma))
    {
        throw std::invalid_argument((boost::format("r0 >= sigma : r0=%.16g, sigma=%.16g") % r0 % sigma).str());
    }

    if (!(t >= 0.0))
    {
        throw std::invalid_argument((boost::format("t >= 0.0 : t=%.16g") % t).str());
    }


    if(t == 0.0)
    {
        return r0;
    }

    const Real psurv(p_survival(t));

    p_int_r_params params = { this, t, rnd * psurv };

    gsl_function F = 
        {
            reinterpret_cast<typeof(F.function)>(&p_int_r_F),
            &params 
        };


    // adjust low and high starting from r0.
    // this is necessary to avoid root finding in the long tails where
    // numerics can be unstable.

    Real low(r0);
    Real high(r0);

    const Real sqrt6Dt(sqrt(6.0 * D * t));
    if(GSL_FN_EVAL(&F, r0) < 0.0)
    {
        // low = r0
        unsigned int H(3);

        for (;;)
        {
            high = r0 + H * sqrt6Dt;

            const Real value(GSL_FN_EVAL(&F, high));
            if(value > 0.0)
            {
                break;
            }

            ++H;

            if(H > 20)
            {
                throw std::runtime_error("drawR: H > 20 while adjusting upper bound of r");
            }
        }

    }
    else
    {
        // high = r0
        unsigned int H(3);

        for (;;)
        {
            low = r0 - H * sqrt6Dt;
            if(low < sigma)
            {
                if(GSL_FN_EVAL(&F, sigma) > 0.0)
                {
                    log_.info("drawR: p_int_r(sigma) > 0.0. "
                              "returning sigma.");
                    return sigma;
                }

                low = sigma;
                break;
            }

            const Real value(GSL_FN_EVAL(&F, low));
            if(value < 0.0)
            {
                break;
            }

            ++H;
        }
    }


    // root finding by iteration.

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
                                                  this->TOLERANCE));

        if(status == GSL_CONTINUE)
        {
            if(i >= maxIter)
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
GreensFunction3DRadInf::Rn(unsigned int n, Real r, Real t,
                            gsl_integration_workspace* workspace,
                            Real tol) const
{
    Real integral;
    Real error;

    p_corr_R_params params = { this, n, r, t };
    gsl_function F = 
        {
            reinterpret_cast<typeof(F.function)>(&p_corr_R_F),
            &params
        };

    const Real umax(sqrt(40.0 / (this->getD() * t))); 

    gsl_integration_qag(&F, 0.0,
                         umax,
                         tol,
                         THETA_TOLERANCE,
                         2000, GSL_INTEG_GAUSS61,
                         workspace, &integral, &error);

    return integral;
}


Real GreensFunction3DRadInf::p_corr_n(unsigned int n, RealVector const& RnTable, RealVector const& lgndTable) const
{
    return RnTable[n] * lgndTable[n] * (2.0 * n + 1.0);
}

Real GreensFunction3DRadInf::ip_corr_n(unsigned int n, RealVector const& RnTable, RealVector const& lgndTable) const
{
    // lgndTable1 is offset by 1; lgndTable1[0] is for n=-1.

    const Real lgnd_n_m1(lgndTable[n]);   // n-1
    const Real lgnd_n_p1(lgndTable[n+2]); // n+1
    
    return RnTable[n] * (lgnd_n_m1 - lgnd_n_p1);// / (1.0 + 2.0 * n);
}


Real GreensFunction3DRadInf::p_corr_table(Real theta, Real r, Real t, RealVector const& RnTable) const
{
    const Index tableSize(RnTable.size());
    if(tableSize == 0)
    {
        return 0.0;
    }

    Real result(0.0);


    Real sin_theta;
    Real cos_theta;
    sincos(theta, &sin_theta, &cos_theta);

    RealVector lgndTable(tableSize);
    gsl_sf_legendre_Pl_array(tableSize-1, cos(theta), &lgndTable[0]);


    const Real p(funcSum_all(boost::bind(&GreensFunction3DRadInf::
                                            p_corr_n,
                                            this,
                                            _1, RnTable, lgndTable),
                               tableSize));

    result = - p * sin_theta;

    result /= 4.0 * M_PI * sqrt(r * r0);

    return result;
}


Real GreensFunction3DRadInf::ip_corr_table(Real theta, Real r,
                                            Real t, RealVector const& RnTable) const
{
    const Index tableSize(RnTable.size());
    if(tableSize == 0)
    {
        return 0.0;
    }

    const Real cos_theta(cos(theta));
    
    // lgndTable is offset by 1. lengTable[0] -> n = -1

    RealVector lgndTable(tableSize + 2);
    lgndTable[0] = 1.0; // n = -1
    gsl_sf_legendre_Pl_array(tableSize, cos_theta, &lgndTable[1]);

    const Real p(funcSum_all(boost::bind(&GreensFunction3DRadInf::
                                            ip_corr_n,
                                            this,
                                            _1, RnTable, lgndTable),
                               tableSize));

    const Real result(- p / (4.0 * M_PI * sqrt(r * r0)));
    return result;
}

Real 
GreensFunction3DRadInf::ip_free(Real theta, Real r, Real t) const
{
    return ip_theta_free(theta, r, r0, t, getD());
}


Real GreensFunction3DRadInf::p_theta(Real theta, Real r, Real t) const
{
    RealVector RnTable;
    makeRnTable(RnTable, r, t);

    return p_theta_table(theta, r, t, RnTable);
}

Real GreensFunction3DRadInf::ip_theta(Real theta, Real r, Real t) const
{
    RealVector RnTable;
    makeRnTable(RnTable, r, t);

    return ip_theta_table(theta, r, t, RnTable);
}


Real GreensFunction3DRadInf::p_theta_table(Real theta, Real r,
                                            Real t, RealVector const& RnTable) const
{
    const Real p_free(this->p_free(theta, r, t));
    const Real p_corr(this->p_corr_table(theta, r, t, RnTable)); 

//    return p_free;
    return (p_free + p_corr);
}

Real GreensFunction3DRadInf::ip_theta_table(Real theta, Real r, Real t, RealVector const& RnTable) const
{
    const Real p_free(this->ip_free(theta, r, t));
    const Real p_corr(this->ip_corr_table(theta, r, t, RnTable)); 

    return (p_free + p_corr);
}

static const Real p_free_max(Real r, Real r0, Real t, Real D)
{
    const Real Dt4(4.0 * D * t);
    const Real Dt4Pi(Dt4 * M_PI);

    const Real term1(exp(- gsl_pow_2(r - r0) / Dt4));
    const Real term2(1.0 / sqrt(Dt4Pi * Dt4Pi * Dt4Pi));

    return term1 * term2;
}

void GreensFunction3DRadInf::makeRnTable(RealVector& RnTable,
                                          Real r, Real t) const
{
    RnTable.clear();

    const Real sigma(getSigma());
    const Real D(getD());
    const Real kf(getkf());

    {  
        // First, estimate the size of p_corr, and if it's small enough,
        // we don't need to calculate it in the first place.
        const Real pirr(p_irr(r, t, r0, kf, D, sigma));
        const Real ipfree_max(ip_free(M_PI, r, t) * 2 * M_PI * r * r);
        
        if(fabs((pirr - ipfree_max) / ipfree_max) < 1e-8)
        {
            return;
        }
    }


    const Real pfreemax(p_free_max(r, r0, t, D));

    gsl_integration_workspace* 
        workspace(gsl_integration_workspace_alloc(2000));
    
    Real Rn_prev(0.0);
    const Real RnFactor(1.0 / (4.0 * M_PI * sqrt(r * r0)));

    const Real integrationTolerance(pfreemax / RnFactor * THETA_TOLERANCE);
    const Real truncationTolerance(pfreemax * THETA_TOLERANCE * 1e-1);
    
    unsigned int n(0);
    for (;;) 
    {
        const Real Rn(this->Rn(n, r, t, workspace, 
                                 integrationTolerance));
        
        RnTable.push_back(Rn);
        
        // truncate when converged enough.
        const Real absRn(fabs(Rn));
        if(absRn * RnFactor < truncationTolerance &&
            absRn < Rn_prev)
        {
            break;
        }

        if(n >= this->MAX_ORDER)
        {
            log_.info("GreensFunction3DRadInf: Rn didn't converge");
            break;
        }
        
        Rn_prev = fabs(Rn);
        
        ++n;
    }

    gsl_integration_workspace_free(workspace);
}

struct GreensFunction3DRadInf::p_theta_params 
{ 
    const GreensFunction3DRadInf* const gf;
    const Real r;
    const Real t;
    GreensFunction3DRadInf::RealVector const& RnTable;
    const Real value;
};
    
Real GreensFunction3DRadInf::ip_theta_F(Real theta, p_theta_params* params)
{
    const GreensFunction3DRadInf* const gf(params->gf); 
    const Real r(params->r);
    const Real t(params->t);
    GreensFunction3DRadInf::RealVector const& RnTable(params->RnTable);
    const Real value(params->value);

    return gf->ip_theta_table(theta, r, t, RnTable) - value;
}


Real GreensFunction3DRadInf::drawTheta(Real rnd, Real r, Real t) const
{
    Real theta;

    const Real sigma(this->getSigma());

    // input parameter range checks.
    if (!(rnd < 1.0 && rnd >= 0.0))
    {
        throw std::invalid_argument((boost::format("rnd < 1.0 && rnd >= 0.0 : rnd=%.16g") % rnd).str());
    }

    if (!(r >= sigma))
    {
        throw std::invalid_argument((boost::format("r >= sigma : r=%.16g, sigma=%.16g") % r % sigma).str());
    }

    if (!(r0 >= sigma))
    {
        throw std::invalid_argument((boost::format("r0 >= sigma : r0=%.16g, sigma=%.16g") % r0 % sigma).str());
    }

    if (!(t >= 0.0))
    {
        throw std::invalid_argument((boost::format("t >= 0.0 : t=%.16g") % t).str());
    }


    // t == 0 means no move.
    if(t == 0.0)
    {
        return 0.0;
    }

    RealVector RnTable;
    makeRnTable(RnTable, r, t);


    // root finding with the integrand form.

    const Real ip_theta_pi(ip_theta_table(M_PI, r, t, RnTable));

    p_theta_params params = { this, r, t, RnTable, rnd * ip_theta_pi };

    gsl_function F = 
        {
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
        const int status(gsl_root_test_interval(low, high, 1e-15,
                                                  THETA_TOLERANCE));

        if(status == GSL_CONTINUE)
        {
            if(i >= maxIter)
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



//
// debug
//

std::string GreensFunction3DRadInf::dump() const
{
    std::ostringstream ss;
    ss << "D = " << this->getD() << ", sigma = " << this->getSigma() <<
        ", kf = " << this->getkf() <<
        ", kD = " << this->getkD() <<
        ", alpha = " << this->getalpha() << std::endl;
    return ss.str();
}

Logger& GreensFunction3DRadInf::log_(
        Logger::get_logger("GreensFunction3DRadInf"));


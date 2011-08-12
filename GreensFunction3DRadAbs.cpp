#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "compat.h"

#include <stdexcept>
#include <vector>
#include <sstream>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_integration.h>

#include "factorial.hpp"
#include "funcSum.hpp"
#include "findRoot.hpp"
#include "freeFunctions.hpp"
#include "SphericalBesselGenerator.hpp"
#include "GreensFunction3DRadAbs.hpp"

const Real GreensFunction3DRadAbs::TOLERANCE;
const Real GreensFunction3DRadAbs::MIN_T_FACTOR;
const unsigned int GreensFunction3DRadAbs::MAX_ORDER;
const unsigned int GreensFunction3DRadAbs::MAX_ALPHA_SEQ;


GreensFunction3DRadAbs::GreensFunction3DRadAbs(
    Real D, Real kf, Real r0, Real Sigma, Real a)
    : GreensFunction3DRadAbsBase(D, kf, r0, Sigma),
      h(kf / (4.0 * M_PI * Sigma * Sigma * D)),
      hsigma_p_1(1.0 + h * Sigma),
      a(a)
{
    const Real sigma(this->getSigma());

    if (a < sigma)
    {
        throw std::invalid_argument((boost::format("a >= sigma : a=%.16g, sigma=%.16g") % a % sigma).str());
    }
    clearAlphaTable();
}

GreensFunction3DRadAbs::~GreensFunction3DRadAbs()
{
    ; // do nothing
}

//
// Alpha-related methods
//

void GreensFunction3DRadAbs::clearAlphaTable() const
{
    std::for_each(this->alphaTable.begin(), this->alphaTable.end(),
                   boost::mem_fn(&RealVector::clear));
    this->alphaOffsetTable[0] = 0;
    std::fill(this->alphaOffsetTable.begin()+1, this->alphaOffsetTable.end(),
               -1);

}


Real GreensFunction3DRadAbs::f_alpha0(Real alpha) const
{
    const Real a(geta());
    const Real sigma(getSigma());

    const Real alpha_a_m_sigma(alpha * (a - sigma));
    const Real hsigma_p_1(this->hsigma_p_1);

    Real sin_alpha_a_m_sigma;
    Real cos_alpha_a_m_sigma;
    sincos(alpha_a_m_sigma, &sin_alpha_a_m_sigma, &cos_alpha_a_m_sigma);

    const Real term1(alpha * sigma * cos_alpha_a_m_sigma);
    const Real term2(hsigma_p_1 * sin_alpha_a_m_sigma);

    const Real result(term1 + term2);

    return result;
}


Real 
GreensFunction3DRadAbs::f_alpha0_aux(Real alpha) const
{
    const Real a(this->geta());
    const Real sigma(this->getSigma());

    const Real term1((a - sigma) * alpha);

    const Real angle(this->hsigma_p_1 / (sigma * alpha));
    const Real term2(std::atan(angle));

    const Real result(term1 - term2);

    return result;
}

struct f_alpha0_aux_params
{ 
    GreensFunction3DRadAbs const* const gf;
    const Real value;
};

static Real f_alpha0_aux_F(Real alpha, f_alpha0_aux_params const* params)
{
    return params->gf->f_alpha0_aux(alpha) - params->value;
}


Real GreensFunction3DRadAbs::alpha0_i(Integer i) const
{
    if (!(i >= 0))
    {
        throw std::out_of_range((boost::format("i >= 0 : i=%.16g") % i).str());
    }


    const Real a(this->geta());
    const Real sigma(this->getSigma());


    const Real target(i * M_PI + M_PI_2);
    f_alpha0_aux_params params = { this, target };


    gsl_function F = 
        { reinterpret_cast<typeof(F.function)>(&f_alpha0_aux_F), &params };


    // We know the range of the solution from - Pi/2 <= atan <= Pi/2.
    const Real interval(M_PI / (a - sigma));
    Real low(i * interval + std::numeric_limits<Real>::epsilon());
    Real high((i+1) * interval);

    //assert(GSL_FN_EVAL(&F, low) * GSL_FN_EVAL(&F, high) < 0.0);

    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));
    gsl_root_fsolver_set(solver, &F, low, high);

    const unsigned int maxIter(100);

    unsigned int j(0);
    for (;;)
    {
        gsl_root_fsolver_iterate(solver);

        low = gsl_root_fsolver_x_lower(solver);
        high = gsl_root_fsolver_x_upper(solver);
        const int status(gsl_root_test_interval(low, high, 0.0, 1e-15));

        if (status == GSL_CONTINUE)
        {
            if (j >= maxIter)
            {
                gsl_root_fsolver_free(solver);
                throw std::runtime_error("alpha0_i: failed to converge");
            }
        }
        else
        {
            break;
        }

        ++j;
    }

    const Real alpha(gsl_root_fsolver_root(solver));
    gsl_root_fsolver_free(solver);
  
    return alpha;
}


void
GreensFunction3DRadAbs::updateAlphaTable0(const Real t) const
{
    RealVector& alphaTable_0(this->getAlphaTable(0));
    alphaTable_0.clear();
    alphaTable_0.reserve(MAX_ALPHA_SEQ);

    const Real alpha0_0(this->alpha0_i(0));
    alphaTable_0.push_back(alpha0_0);

    const Real Dt(this->getD() * t);

//    const Real alpha_cutoff(sqrt((- log(TOLERANCE * 1e-2) / Dt)
//                                 + alpha0_0 * alpha0_0));
    const Real alpha_cutoff(sqrt((- log(TOLERANCE * 1e-3) / Dt)));

    unsigned int i(1);
    for (;;)
    {
        const Real alpha0_i(this->alpha0_i(i));
        alphaTable_0.push_back(alpha0_i);

        if (alpha0_i > alpha_cutoff && i >= 10) // make at least 10 terms
        {
            break;
        }

        ++i;

        if (i >= MAX_ALPHA_SEQ)
        {
            break;
        }
    }
}

Real GreensFunction3DRadAbs::f_alpha(Real alpha, Integer n) const
{
    const Real a(this->geta());
    const Real sigma(getSigma());
    const Real aAlpha(a * alpha);
    const Real sigmaAlpha(getSigma() * alpha);
    const Real hSigma(geth() * getSigma());
    const Real realn(static_cast<Real>(n));

    const Real hSigma_m_n(hSigma - realn);


    const SphericalBesselGenerator& s(SphericalBesselGenerator::instance());

    const Real js1(s.j(n,   sigmaAlpha));
    const Real ys1(s.y(n,   sigmaAlpha));
    const Real js2(s.j(n+1, sigmaAlpha));
    const Real ys2(s.y(n+1, sigmaAlpha));
    const Real ja( s.j(n,   aAlpha));
    const Real ya( s.y(n,   aAlpha));

    const Real term1((hSigma_m_n * js1 + sigmaAlpha * js2) * ya);
    const Real term2((hSigma_m_n * ys1 + sigmaAlpha * ys2) * ja);

    const Real factor(2.0 * alpha * sqrt(a * sigma) * M_1_PI);

    const Real result((term1 - term2) * factor);
    
    return result;
}

static inline const Real G(const unsigned int n, const unsigned int k)
{
    return factorial(n + k) * (factorial_r(k) * factorial_r(n - k));
}


static Real P(Integer n, Real x)
{
    Real result(0.0);

    Real sx2(1.0);
    Integer term1(1);

    const Real x2sq_r(1.0 / gsl_pow_2(x + x));
    const unsigned int maxm(n / 2);
    for(unsigned int m(0); m <= maxm; ++m)
    {
        const Real value(term1 * sx2 * G(n, 2 * m));
        result += value;

        term1 = - term1;
        sx2 *= x2sq_r;
    }

    return result;
}

static std::pair<Real, Real> P2(Integer n, Real x)
{
    Real result(0.0);
    Real resultp(0.0);

    Real sx2(1.0);
    Integer term1(1);

    const Real x2sq_r(1.0 / gsl_pow_2(x + x));
    const unsigned int np1(n + 1);
    const unsigned int maxm(n / 2);
    for(unsigned int m(0); m <= maxm; ++m)
    {
        const Real sx2p(term1 * sx2);
        const unsigned int m2(2 * m);
        const Real value(sx2p * G(n, m2));
        result += value;

        const Real valuep(sx2p * G(np1, m2));
        resultp += valuep;

        term1 = - term1;
        sx2 *= x2sq_r;
    }

    if (n % 2)
    {
        resultp += term1 * sx2 * G(np1, np1);
    }


    return std::make_pair(result, resultp);
}


static Real Q(Integer n, Real x)
{
    Real result(0.0);

    Real sx2(1.0 / (x + x));
    Integer term1(1);

    const Real x2sq(sx2 * sx2);
    const unsigned int maxm((n+1)/2); // sum_(0)^((n-1)/2)
    for(unsigned int m(0); m < maxm; ++m)
    {
        const Real value(term1 * sx2 * G(n, 2 * m + 1));
        result += value;

        term1 = - term1;  // (-1)^m
        sx2 *= x2sq;
    }

    return result;
}

static std::pair<Real, Real> Q2(Integer n, Real x)
{
    Real result(0.0);
    Real resultp(0.0);

    Real sx2(1.0 / (x + x));
    Integer term1(1);  // (-1)^m

    const Real x2sq(sx2 * sx2);
    const unsigned int np1(n + 1);
    const unsigned int maxm((n+1)/2); // sum_(0)^((n-1)/2)
    for(unsigned int m(0); m < maxm; ++m)
    {
        const Real sx2p(term1 * sx2);
        const unsigned int m2p1(2 * m + 1);
        const Real value(sx2p * G(n, m2p1));
        result += value;

        const Real valuep(sx2p * G(np1, m2p1));
        resultp += valuep;

        term1 = - term1; // (-1)^m
        sx2 *= x2sq;
    } 


    if (!(n % 2))
    {
        resultp += term1 * sx2 * G(np1, np1);
    }


    return std::make_pair(result, resultp);
}


Real GreensFunction3DRadAbs::f_alpha_aux(Real alpha, Integer n) const
{
    if (alpha == 0.0)
    {
        return -1.0;
    }

    const Real a(geta());
    const Real sigma(getSigma());

    const Real aAlpha(a * alpha);
    const Real sigmaAlpha(sigma * alpha);

    const Real n_m_hSigma(n - h * sigma);

    /*(a - s) u - 
      ArcTan[(P[n, a u] ((-n + h s) P[n, s u] + s u Q[1 + n, s u]) -
               Q[n, a u] (s u P[1 + n, s u] + (n - h s) Q[n, s u]))/
             (Q[n, a u] ((-n + h s) P[n, s u] + s u Q[1 + n, s u]) + 
               P[n, a u] (s u P[1 + n, s u] + (n - h s) Q[n, s u]))]
    */

    const Real Pa(P(n, aAlpha));
    const Real Qa(Q(n, aAlpha));

    Real Ps;
    Real Psp;
    // boost::tie(Ps, Psp) = P2(n, sigmaAlpha);
    {
        std::pair<Real,Real> res(P2(n, sigmaAlpha));
        Ps = res.first;
        Psp = res.second;
    }

    Real Qs;
    Real Qsp;
    // boost::tie(Qs, Qsp) = Q2(n, sigmaAlpha);
    {
        std::pair<Real,Real> res(Q2(n, sigmaAlpha));
        Qs = res.first;
        Qsp = res.second;
    }

    const Real n_m_hSigmaPs(n_m_hSigma * Ps);
    const Real n_m_hSigmaQs(n_m_hSigma * Qs);
    const Real sigmaAlphaPsp(sigmaAlpha * Psp);
    const Real sigmaAlphaQsp(sigmaAlpha * Qsp);

    const Real Qa_Pa(Qa / Pa);

    const Real A(sigmaAlphaQsp - n_m_hSigmaPs);
    const Real B(sigmaAlphaPsp + n_m_hSigmaQs);

    // this form, dividing all terms by Pa, prevents overflow.
    const Real angle((A - Qa_Pa * B) / (Qa_Pa * A + B));

    const Real term1((a - sigma) * alpha);
    const Real term2(std::atan(angle));

    const Real result(term1 - term2);

    return result;
}

struct f_alpha_aux_params
{ 
    GreensFunction3DRadAbs const* const gf;
    const Integer n;
    const Real value;
};

static Real f_alpha_aux_F(Real alpha, f_alpha_aux_params const* params)
{
    return params->gf->f_alpha_aux(alpha, params->n) - params->value;
}

Real 
GreensFunction3DRadAbs::alpha_i(Integer i, Integer n, 
                                        gsl_root_fsolver* solver) const
{
    const Real sigma(this->getSigma());
    const Real a(this->geta());

    const Real target(M_PI * i + M_PI_2);

    const Real factor(1.0 / (a - sigma));
    Real low((target - M_PI_2) * factor);
    Real high((target + M_PI_2) * factor);

    f_alpha_aux_params params = { this, n, target };

    gsl_function F = 
        { reinterpret_cast<typeof(F.function)>(&f_alpha_aux_F), &params };

    gsl_root_fsolver_set(solver, &F, low, high);

    const unsigned int maxIter(100);
    unsigned int k(0);
    for (;;)
    {
        gsl_root_fsolver_iterate(solver);
        
        low = gsl_root_fsolver_x_lower(solver);
        high = gsl_root_fsolver_x_upper(solver);
        const int status(gsl_root_test_interval(low, high, 1e-6, 1e-15));
        
        if (status == GSL_CONTINUE)
        {
            if (k >= maxIter)
            {
                gsl_root_fsolver_free(solver);
                throw std::runtime_error("alpha_i: failed to converge");
            }
        }
        else
        {
            break;
        }
        
        ++k;
    }
    
    const Real alpha(gsl_root_fsolver_root(solver));

    return alpha;
}


unsigned int
GreensFunction3DRadAbs::alphaOffset(unsigned int n) const
{
    if (this->alphaOffsetTable[n] >= 0)
    {
        return this->alphaOffsetTable[n];
    }

    const Real sigma(this->getSigma());
    const Real a(this->geta());

    assert(this->alphaOffsetTable.size() >= n);
    unsigned int offset(this->alphaOffsetTable[n-1]);

    const Real factor(1.0 / (a - sigma));

    Real target(offset * M_PI + M_PI_2);
    // We know the range of the solution from - Pi/2 <= atan <= Pi/2.
    const Real alphaMid(target * factor);
    const Real alphaHalfRange(M_PI_2 * factor);
    Real low(alphaMid - alphaHalfRange * (1.0 - 1e-3)); // avoid zero.
    Real high(alphaMid + alphaHalfRange);


    // Here we find the interval where the first positive root is in.
    // We find the first pair of alpha
    // (Pi * offset + Pi/2) +- Pi/2 / (a - sigma)
    // where the values of f_alpha() straddle.
    // The assumption is the interval between roots is not much
    // smaller than Pi / (a - sigma).


    Real lowvalue(f_alpha(low,n));
    Real highvalue(f_alpha(high,n));

    for (;;) // this can be much faster if better initial guess is given.
    {

        if (lowvalue * highvalue < 0) // low and high straddle?
        {
            break;
        }

        ++offset;
        target = M_PI * offset + M_PI_2;
        low = (target - M_PI_2) * factor;
        high = (target + M_PI_2) * factor;

        lowvalue = highvalue;
        highvalue = f_alpha(high, n);
    }

    this->alphaOffsetTable[n] = offset;

    return offset;
}


void
GreensFunction3DRadAbs::updateAlphaTable(const unsigned int n,
                                                  const Real t) const
{
    if (!(n >= 0 && n <= this->MAX_ORDER))
    {
        throw std::range_error((boost::format("n >= 0 && n <= this->MAX_ORDER : n=%.16g, this->MAX_ORDER=%.16g") % n % this->MAX_ORDER).str());
    }


    if (n == 0)
    {
        this->updateAlphaTable0(t);
        return;
    }

    const unsigned int offset(alphaOffset(n));

    RealVector& alphaTable_n(this->getAlphaTable(n));
    alphaTable_n.clear();
    alphaTable_n.reserve(MAX_ALPHA_SEQ);

    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));

    const Real alphan_0(alpha_i(offset, n, solver));
    const Real alphan_0_sq(alphan_0 * alphan_0);

    alphaTable_n.push_back(alphan_0);

    const Real Dt(this->getD() * t);

    const Real threshold(this->TOLERANCE * 1e-2 * 
                          alphan_0_sq * exp(- Dt * alphan_0_sq));
   
    const unsigned int end(offset + MAX_ALPHA_SEQ);
    unsigned int i(offset + 1);
    for (;;)
    {
        const Real alpha_i(this->alpha_i(i, n, solver));

        alphaTable_n.push_back(alpha_i);

        // cutoff
        const Real alpha_i_sq(alpha_i * alpha_i);
        if (alpha_i_sq * exp(- Dt * alpha_i_sq)  < threshold)
        {
            break;
        }


        ++i;

        if (i >= end)
        {
            log_.info("alphaTable (%d): didn't converge. t = %.16g, %s",
                       n, t, dump().c_str());
            break;
        }
    }

    gsl_root_fsolver_free(solver);
}




Real 
GreensFunction3DRadAbs::p_0_i(Real alpha, Real r) const
{
    const Real a(geta());
    const Real sigma(getSigma());
    const Real h(geth());
    const Real hsigma_p_1(this->hsigma_p_1);

    const Real sigmasq(sigma * sigma);
    const Real alphasq(alpha * alpha);

    Real num1;
    {
        const Real angle_r(alpha * (r - sigma));
        Real sin_r;
        Real cos_r;
        sincos(angle_r, &sin_r, &cos_r);
        num1 = alpha * sigma * cos_r + hsigma_p_1 * sin_r ;
    }

    const Real num2(num_r0(alpha));

    const Real den(2 * M_PI * r * r0 * 
                    ((a - sigma) * sigmasq * alphasq +
                      hsigma_p_1 * (a + a * h * sigma - h * sigmasq)));

    const Real result(num1 * num2 / den);

    return result;
}


Real 
GreensFunction3DRadAbs::p_survival_i(Real alpha) const
{
    const Real a(geta());
    const Real sigma(getSigma());
    const Real h(geth());
    const Real hsigma_p_1(this->hsigma_p_1);

    const Real sigmasq(sigma * sigma);
    const Real alphasq(alpha * alpha);

    const Real angle_a(alpha * (a - sigma));
    const Real cos_a(cos(angle_a));

    const Real num1(h * sigmasq * hsigma_p_1 
                     - a * (hsigma_p_1 * hsigma_p_1
                             + sigmasq * alphasq) * cos_a);

    const Real num2(num_r0(alpha));

    const Real den(r0 * hsigma_p_1 * alpha * 
                    (- hsigma_p_1 *
                      (a + a * h * sigma - h * sigmasq) 
                      + (sigma - a) * sigmasq * alphasq));

    const Real result(- 2.0 * num1 * num2 / den);

    return result;
}


Real 
GreensFunction3DRadAbs::dp_survival_i(Real alpha) const
{
    const Real a(geta());
    const Real sigma(getSigma());
    const Real h(geth());
    const Real hsigma_p_1(this->hsigma_p_1);

    const Real sigmasq(sigma * sigma);
    const Real alphasq(alpha * alpha);

    const Real angle_a(alpha * (a - sigma));
    const Real cos_a(cos(angle_a));

    const Real num1(alpha * (h * sigmasq * hsigma_p_1 
                               - (a * (hsigma_p_1 * hsigma_p_1 
                                         + sigmasq * alphasq)) * cos_a));

    const Real num2(num_r0(alpha));

    const Real den(r0 * hsigma_p_1 * 
                    (- hsigma_p_1 * (a + a * h * sigma - h * sigmasq))
                    + (sigma - a) * sigmasq * alphasq);

    const Real result(2.0 * getD() * num1 * num2 / den);

    return result;
}


Real 
GreensFunction3DRadAbs::leavea_i(Real alpha) const
{
    const Real a(geta());
    const Real sigma(getSigma());
    const Real h(geth());
    const Real D(getD());
    const Real hsigma_p_1(this->hsigma_p_1);

    const Real sigmasq(sigma * sigma);
    const Real alphasq(alpha * alpha);

    const Real angle_a(alpha * (a - sigma));
    const Real cos_a(cos(angle_a));

    const Real num1(alpha * (hsigma_p_1 * hsigma_p_1 + sigmasq * alphasq)
                     * cos_a);

    const Real num2(num_r0(alpha));
    
    const Real den(2 * a * M_PI * r0 * hsigma_p_1 *
                    (hsigma_p_1 * (a + a * h * sigma - h * sigmasq)
                      + (a - sigma) * sigmasq * alphasq));

    const Real result(D * num1 * num2 / den);

    return result;
}

Real GreensFunction3DRadAbs::leaves_i(Real alpha) const
{
    const Real a(geta());
    const Real sigma(getSigma());
    const Real h(geth());
    const Real D(getD());
    const Real hsigma_p_1(this->hsigma_p_1);

    const Real sigmasq(sigma * sigma);
    const Real alphasq(alpha * alpha);

    const Real num(h * alpha * num_r0(alpha));
                      
    const Real den(2 * M_PI * r0 *
                    ((a - sigma) * sigmasq * alphasq +
                      hsigma_p_1 * (a + a * h * sigma - h * sigmasq)));

    const Real result(- D * num / den);
        
    return result;
}


Real GreensFunction3DRadAbs::p_leavea_i(Real alpha,
                                                Real pleave_factor) const
{
    const Real a(geta());
    const Real sigma(getSigma());

    const Real hsigma_p_1(this->hsigma_p_1);
    const Real sigmasq(sigma * sigma);
    const Real alphasq(alpha * alpha);

    const Real angle_a(alpha * (a - sigma));
    const Real cos_a(cos(angle_a));

    const Real num1((hsigma_p_1 * hsigma_p_1 + sigmasq * alphasq) * cos_a);

    const Real result(- 2.0 * a * num1 * pleave_factor / hsigma_p_1);

    return result;
}


Real 
GreensFunction3DRadAbs::p_leaves_i(Real alpha,
                                           Real pleave_factor) const
{
    const Real sigma(getSigma());
    const Real h(geth());
 
    const Real num(h * sigma * sigma);
                      
    const Real result(2.0 * num * pleave_factor);
        
    return result;
}

Real 
GreensFunction3DRadAbs::p_survival_den(Real alpha) const
{
    const Real a(geta());
    const Real sigma(getSigma());
    const Real h(geth());
    const Real hsigma_p_1(this->hsigma_p_1);
    const Real sigmasq(sigma * sigma);
    const Real alphasq(alpha * alpha);

    const Real den(r0 * alpha *
                    ((a - sigma) * sigmasq * alphasq +
                      hsigma_p_1 * (a + a * h * sigma - h * sigmasq)));
    
    return den;
}



Real GreensFunction3DRadAbs::num_r0(Real alpha) const
{
    const Real sigma(getSigma());
    const Real angle_r0(alpha * (r0 - sigma));

    Real sin_r0;
    Real cos_r0;
    sincos(angle_r0, &sin_r0, &cos_r0);

    const Real hsigma_p_1(this->hsigma_p_1);
    const Real result(alpha * sigma * cos_r0 + hsigma_p_1 * sin_r0);

    return result;
}


Real GreensFunction3DRadAbs::pleaveFactor(Real alpha) const
{
    return num_r0(alpha) / p_survival_den(alpha);
}


Real
GreensFunction3DRadAbs::p_int_r_i(Real r, Real alpha,
                                          Real num_r0) const
{
    const Real sigma(getSigma());

    const Real angle_r(alpha * (r - sigma));
    Real sin_r;
    Real cos_r;
    sincos(angle_r, &sin_r, &cos_r);  // do sincos here; latency. 

    const Real h(geth());
    const Real hsigma_p_1(this->hsigma_p_1);

    const Real sigmasq(sigma * sigma);
    const Real alphasq(alpha * alpha);

    const Real hsigma(h * sigma);

    const Real num1(alpha * (hsigma * sigma - hsigma * r * cos_r 
                               - (r - sigma) * cos_r) 
                     + (hsigma_p_1 + r * sigma * alphasq) * sin_r);

    const Real num2(num_r0);

    const Real den(r0 * alphasq * 
                    ((a - sigma) * sigmasq * alphasq +
                      hsigma_p_1 * (a + a * h * sigma - h * sigmasq)));

    const Real result(2 * num1 * num2 / den);

    return result;
}


void 
GreensFunction3DRadAbs::createPsurvTable(RealVector& table) const
{
    const RealVector& alphaTable_0(this->getAlphaTable(0));

    table.clear();
    table.reserve(alphaTable_0.size());

    std::transform(alphaTable_0.begin(), alphaTable_0.end(),
                    std::back_inserter(table),
                    boost::bind(&GreensFunction3DRadAbs::p_survival_i,
                                 this, _1));
}


void 
GreensFunction3DRadAbs::createNum_r0Table(RealVector& table) const
{
    const RealVector& alphaTable_0(this->alphaTable[0]);

    table.clear();
    table.reserve(alphaTable_0.size());

    std::transform(alphaTable_0.begin(), alphaTable_0.end(),
                    std::back_inserter(table),
                    boost::bind(&GreensFunction3DRadAbs::num_r0,
                                 this, _1));
}

void 
GreensFunction3DRadAbs::createPleaveFactorTable(RealVector& table) const
{
    const RealVector& alphaTable_0(this->alphaTable[0]);

    table.clear();
    table.reserve(alphaTable_0.size());

    std::transform(alphaTable_0.begin(), alphaTable_0.end(),
                    std::back_inserter(table),
                    boost::bind(&GreensFunction3DRadAbs::pleaveFactor,
                                 this, _1));
}


void 
GreensFunction3DRadAbs::createPleavesTable(RealVector& table,
                                                   RealVector const& pleaveFactorTable) const
{
    const RealVector& alphaTable_0(this->alphaTable[0]);

    assert(pleaveFactorTable.size() >= alphaTable_0.size());

    table.clear();
    table.reserve(alphaTable_0.size());

    for(unsigned int i(0); i < alphaTable_0.size(); ++i)
    {
        const Real alpha(alphaTable_0[i]);
        table.push_back(p_leaves_i(alpha, pleaveFactorTable[i]));
    }
}

void 
GreensFunction3DRadAbs::createPleaveaTable(RealVector& table,
                                                   RealVector const& pleaveFactorTable) const
{
    const RealVector& alphaTable_0(this->alphaTable[0]);

    assert(pleaveFactorTable.size() >= alphaTable_0.size());

    table.clear();
    table.reserve(alphaTable_0.size());

    for(unsigned int i(0); i < alphaTable_0.size(); ++i)
    {
        const Real alpha(alphaTable_0[i]);
        table.push_back(p_leavea_i(alpha, pleaveFactorTable[i]));
    }
}


Real 
GreensFunction3DRadAbs::p_0_i_exp(unsigned int i, Real t, Real r) const
{
    const Real alpha(this->getAlpha0(i));
    return std::exp(- getD() * t * alpha * alpha) * p_0_i(alpha, r);
}


Real 
GreensFunction3DRadAbs::p_survival_i_exp(unsigned int i, Real t) const
{
    const Real alpha(this->getAlpha0(i));
    return p_survival_i_alpha(alpha, t);
}

Real 
GreensFunction3DRadAbs::p_survival_i_alpha(Real alpha, Real t) const
{
    return std::exp(- getD() * t * alpha * alpha) * 
        p_survival_i(alpha);
}

Real 
GreensFunction3DRadAbs::p_survival_2i_exp(unsigned int i,
                                                  Real t) const
{
    const Real Dt(getD() * t);
    const Real alpha0(this->getAlpha0(2 * i));
    const Real p0(std::exp(- Dt * alpha0 * alpha0) * 
                   p_survival_i(alpha0));

    const Real alpha1(this->getAlpha0(2 * i + 1));
    const Real p1(std::exp(- Dt * alpha1 * alpha1) * 
                   p_survival_i(alpha1));

    return p0 + p1;
}

Real 
GreensFunction3DRadAbs::p_survival_i_exp_table(unsigned int i,
                                                       Real t,
                                                       RealVector const& table) const
{
    const Real alpha(this->getAlpha0(i));
    return std::exp(- getD() * t * alpha * alpha) * table[i];
}

Real 
GreensFunction3DRadAbs::p_leave_i_exp_table(unsigned int i, Real t, RealVector const& table) const
{
    const Real alpha(getAlpha0(i));
    return expm1(- getD() * t * alpha * alpha) * table[i];
}


Real 
GreensFunction3DRadAbs::dp_survival_i_exp(unsigned int i, Real t) const
{
    const Real alpha(this->getAlpha0(i));
    return std::exp(- getD() * t * alpha * alpha) * 
        dp_survival_i(alpha);
}

Real 
GreensFunction3DRadAbs::leavea_i_exp(unsigned int i, Real t) const
{
    const Real alpha(this->getAlpha0(i));
    return std::exp(- getD() * t * alpha * alpha) * leavea_i(alpha);
}

Real 
GreensFunction3DRadAbs::leaves_i_exp(unsigned int i, Real t) const
{
    const Real alpha(this->getAlpha0(i));

    return std::exp(- getD() * t * alpha * alpha) * leaves_i(alpha);
}

Real 
GreensFunction3DRadAbs::p_leavea_i_exp(unsigned int i,
                                               Real t) const
{
    const Real alpha(this->getAlpha0(i));
    const Real num_r0(this->num_r0(alpha)); 
    const Real den(this->p_survival_den(alpha)); 
    return exp(- getD() * t * alpha * alpha) * 
        p_leavea_i(alpha, num_r0 / den);
}

Real 
GreensFunction3DRadAbs::p_leaves_i_exp(unsigned int i, Real t) const
{
    const Real alpha(this->getAlpha0(i));
    const Real num_r0(this->num_r0(alpha)); 
    const Real den(this->p_survival_den(alpha)); 
    return exp(- getD() * t * alpha * alpha) * 
        p_leaves_i(alpha, num_r0 / den);
}

Real 
GreensFunction3DRadAbs::p_int_r_i_exp(unsigned int i,
                                              Real t,
                                              Real r) const
{
    const Real alpha(this->getAlpha0(i));

    return std::exp(- getD() * t * alpha * alpha) * 
        p_int_r_i(r, alpha, num_r0(alpha));
}

Real 
GreensFunction3DRadAbs::p_int_r_i_exp_table(unsigned int i,
                                                    Real t,
                                                    Real r,
                                                    RealVector& num_r0Table) const
{
    const Real alpha(this->getAlpha0(i));
    return std::exp(- getD() * t * alpha * alpha) * 
        p_int_r_i(r, alpha, num_r0(alpha));//num_r0Table[i]);
}

Real 
GreensFunction3DRadAbs::p_0(Real t, Real r) const
{
    const Real p(funcSum(boost::bind(&GreensFunction3DRadAbs::
                                        p_0_i_exp,
                                        this,
                                        _1, t, r),
                           this->MAX_ALPHA_SEQ));
    return p;
}


unsigned int
GreensFunction3DRadAbs::guess_maxi(Real t) const
{
    const unsigned int safety(2);

    if (t >= INFINITY)
    {
        return safety;
    }

    const Real D(this->getD());
    const Real sigma(this->getSigma());
    const Real a(this->geta());

    const Real alpha0(this->getAlpha0(0));
    const Real Dt(D * t);

    const Real thr(exp(- Dt * alpha0 * alpha0) * this->TOLERANCE * 1e-1);

    if (thr <= 0.0)
    {
        return this->MAX_ALPHA_SEQ;
    }

    const Real max_alpha(sqrt(alpha0 * alpha0 - log(thr) / Dt));

    const unsigned int 
        maxi(safety + 
              static_cast<unsigned int>(max_alpha * (a - sigma) / M_PI));

    return std::min(maxi, this->MAX_ALPHA_SEQ);
}


Real GreensFunction3DRadAbs::p_survival(Real t) const
{
    RealVector psurvTable;

    const Real p(p_survival_table(t, psurvTable));

    return p;
}

Real 
GreensFunction3DRadAbs::p_survival_table(Real t, RealVector& psurvTable) const
{
    Real p;

    const Real D(this->getD());
    const Real sigma(getSigma());
    const Real a(this->geta());

    const Real distToa(a - r0);
    const Real distTos(r0 - sigma);

    const Real H(6.0); // a fairly strict criterion for safety.
    const Real maxDist(H * sqrt(6.0 * D * t));

    if (distToa > maxDist)
    {
        if (distTos > maxDist) // far from anything; it'll survive.
        {
            p = 1.0;  
        }
        else // close only to s, ignore a
        {
            const Real sigma(this->getSigma());
            const Real kf(this->getkf());
            p = p_survival_irr(t, r0, kf, D, sigma);
        }
    }
    else
    {
        if (distTos > maxDist)  // close only to a.
        {
            p = p_survival_nocollision(t, r0, D, a);
        }
        else  // close to both boundaries.  do the normal calculation.
        {
            const unsigned int maxi(guess_maxi(t));
            
            if (psurvTable.size() < maxi + 1)
            {
                IGNORE_RETURN getAlpha0(maxi);  // this updates the table
                this->createPsurvTable(psurvTable);
            }

            p = funcSum_all(boost::bind(&GreensFunction3DRadAbs::
                                          p_survival_i_exp_table, 
                                          this,
                                          _1, t, psurvTable),
                             maxi);
        }
    }

    return p;
}

Real 
GreensFunction3DRadAbs::p_leave_table(Real t, RealVector const& table) const
{
    return funcSum(
        boost::bind(&GreensFunction3DRadAbs::p_leave_i_exp_table, 
                    this, _1, t, table),
        table.size());
}


Real GreensFunction3DRadAbs::dp_survival(Real t) const
{
    return funcSum(
        boost::bind(&GreensFunction3DRadAbs::dp_survival_i_exp, 
                    this, _1, t),
        MAX_ALPHA_SEQ);
}


Real GreensFunction3DRadAbs::leaves(Real t) const
{
    return funcSum(
        boost::bind(&GreensFunction3DRadAbs::leaves_i_exp,
                    this, _1, t),
        MAX_ALPHA_SEQ);
}

Real GreensFunction3DRadAbs::leavea(Real t) const
{
    return funcSum(
        boost::bind(&GreensFunction3DRadAbs::leavea_i_exp,
                    this, _1, t),
        MAX_ALPHA_SEQ);
}

Real GreensFunction3DRadAbs::p_leaves(Real t) const
{
    return funcSum_all(
        boost::bind(&GreensFunction3DRadAbs::p_leaves_i_exp,
                    this, _1, t),
        guess_maxi(t));
}

Real GreensFunction3DRadAbs::p_leavea(Real t) const
{
    return funcSum_all(
        boost::bind(&GreensFunction3DRadAbs::p_leavea_i_exp,
                    this, _1, t),
        guess_maxi(t));
}

Real GreensFunction3DRadAbs::p_int_r(Real r, Real t) const
{
    return funcSum(
        boost::bind(&GreensFunction3DRadAbs::p_int_r_i_exp,
                    this, _1, t, r),
        MAX_ALPHA_SEQ);
}

Real GreensFunction3DRadAbs::p_int_r_table(Real r, Real t, RealVector const& num_r0Table) const
{
    return funcSum(
        boost::bind(&GreensFunction3DRadAbs::p_int_r_i_exp_table,
                    this, _1, t, r, num_r0Table), num_r0Table.size());
}

struct p_survival_table_params
{ 
    GreensFunction3DRadAbs const* const gf;
    GreensFunction3DRadAbs::RealVector& table;
    const Real rnd;
};

Real p_survival_table_F(Real t, p_survival_table_params const* params)
{
    return params->rnd - params->gf->p_survival_table(t, params->table);
}

struct p_survival_params
{ 
    GreensFunction3DRadAbs const* const gf;
    const Real rnd;
};

static Real p_survival_F(Real t, p_survival_params const* params)
{
    return params->rnd - params->gf->p_survival(t);
}

struct p_survival_2i_params
{ 
    GreensFunction3DRadAbs const* const gf;
    const Real t;
};

static Real p_survival_2i_F(Real ri, p_survival_2i_params const* params)
{
    return params->gf->p_survival_2i_exp(static_cast<unsigned int>(ri),
                                         params->t);
}

struct p_survival_i_alpha_params
{ 
    GreensFunction3DRadAbs const* const gf;
    const Real t;
};

static Real p_survival_i_alpha_F(Real alpha,
                                 p_survival_i_alpha_params const* params)
{
    return params->gf->p_survival_i_alpha(alpha, params->t);
}

struct p_leave_params
{ 
    GreensFunction3DRadAbs const* const gf;
    GreensFunction3DRadAbs::RealVector const& table;
    const Real rnd;
};

Real p_leave_F(Real t, p_leave_params const* params)
{
    return - params->gf->p_leave_table(t, params->table) - params->rnd;
}

struct p_int_r_params
{ 
    GreensFunction3DRadAbs const* const gf;
    const Real t;
    const Real rnd;
};


static Real p_int_r_F(Real r, p_int_r_params const* params)
{
    return params->gf->p_int_r(r, params->t) - params->rnd;
}

Real GreensFunction3DRadAbs::drawTime(Real rnd) const
{
    const Real D(this->getD());
    const Real sigma(this->getSigma());
    const Real kf(this->getkf());
    const Real a(this->geta());

    if (!(rnd < 1.0 && rnd >= 0.0))
    {
        throw std::invalid_argument((boost::format("rnd < 1.0 && rnd >= 0.0 : rnd=%.16g") % rnd).str());
    }

    if (!(r0 >= sigma && r0 <= a))
    {
        throw std::invalid_argument((boost::format("r0 >= sigma && r0 <= a : r0=%.16g, sigma=%.16g, a=%.16g") % r0 % sigma % a).str());
    }


    if (r0 == a || a == sigma)
    {
        return 0.0;
    }

    Real t_guess;
    Real dist;

    if (kf != 0)
    {
        dist = std::min(a - r0, r0 - sigma);
    }
    else
    {
        dist = a - r0;
    }

    t_guess = dist * dist / (6.0 * D);
    t_guess *= .1;

    const Real minT(std::min(sigma * sigma / D * this->MIN_T_FACTOR,
                               t_guess * 1e-6));

    RealVector psurvTable;

    p_survival_table_params params = { this, psurvTable, rnd };

    gsl_function F = 
        {
            reinterpret_cast<typeof(F.function)>(&p_survival_table_F),
            &params 
        };

    Real low(t_guess);
    Real high(t_guess);

    // adjust high and low to make sure that f(low) and f(high) straddle.
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

            if (fabs(high) >= 1e10)
            {
                throw std::runtime_error(
                    (boost::format(
                        "couldn't adjust high. F(%.16g) = %.16g; r0 = %.16g, %s") %
                        high % GSL_FN_EVAL(&F, high) % r0 % dump()).str());
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
            
            // FIXME: 
            if (fabs(low) <= minT ||
                fabs(low_value - low_value_prev) < TOLERANCE) 
            {
                log_.info("couldn't adjust low. F(%.16g) = %.16g; r0 = %.16g, %s",
                          low, GSL_FN_EVAL(&F, low), r0,
                          dump().c_str());
                log_.info("returning %.16g", low);
                return low;
            }
            low_value_prev = low_value;

            low *= .1;
        }
    }

    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));

    const Real t(findRoot(F, solver, low, high, 0.0, 
                            TOLERANCE, "drawTime"));

    gsl_root_fsolver_free(solver);

    return t;
}

GreensFunction3DRadAbs::EventKind
GreensFunction3DRadAbs::drawEventType(Real rnd, Real t) const
{
    const Real D(this->getD());
    const Real sigma(this->getSigma());
    const Real kf(this->getkf());
    const Real a(this->geta());

    if (!(rnd < 1.0 && rnd >= 0.0))
    {
        throw std::invalid_argument((boost::format("rnd < 1.0 && rnd >= 0.0 : rnd=%.16g") % rnd).str());
    }

    if (!(r0 >= sigma && r0 < a))
    {
        throw std::invalid_argument((boost::format("r0 >= sigma && r0 < a : r0=%.16g, sigma=%.16g, a=%.16g") % r0 % sigma % a).str());
    }

    if (!(t > 0.0))
    {
        throw std::invalid_argument((boost::format("t > 0.0 : t=%.16g") % t).str());
    }


    if (kf == 0)
    {
        return IV_ESCAPE;
    }
    
    // First, check if r0 is close only either to a or sigma relative
    // to Dt.  In such cases, the event type is always IV_ESCAPE or 
    // IV_REACTION, respectively. This avoids numerical instability in 
    // calculating leavea() and/or leaves().

    // Here, use a rather large threshold for safety.
    const unsigned int H(6); 
    const Real max_dist(H * sqrt(6.0 * D * t));
    const Real a_dist(a - r0);
    const Real s_dist(r0 - sigma);


    if (a_dist > max_dist)
    {
        if (s_dist < max_dist)
        {
            return IV_REACTION;
        }
    }
    else // a_dist < max_dist
    {
        if (s_dist > max_dist)
        {
            return IV_ESCAPE;
        }
    }

    const Real reaction(leaves(t) * 4.0 * M_PI * sigma * sigma);
    const Real escape(leavea(t) * 4.0 * M_PI * a * a);
    const Real value(reaction / (reaction + escape));

    if (rnd <= value)  
    {
        return IV_REACTION;   // leaves
    }
    else 
    {
        return IV_ESCAPE;     // leavea
    }
}

Real 
GreensFunction3DRadAbs::drawPleavea(gsl_function const& F,
                                            gsl_root_fsolver* solver,
                                            Real t_guess,
                                            RealVector& pleaveFactorTable,
                                            RealVector& pleaveaTable) const
{
    const Real minT(1e-12);

    Real low(t_guess);
    Real high(t_guess);

    // adjust high and low to make sure that f(low) and f(high) straddle.
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

            if (fabs(high) >= 1e10)
            {
                throw std::runtime_error(
                    (boost::format(
                        "couldn't adjust high. Fa(%.16g) = %.16g; r0 = %.16g, %s") %
                        high % GSL_FN_EVAL(&F, high) % r0 % dump()).str());
            }

            log_.info("drawTime2: adjusting high: %.16g Fa = %.16g", high, high_value);
            high *= 10;
        }
    }
    else
    {
        Real low_value_prev(value);
        low *= .1;

        for (;;)
        {
            this->updateAlphaTable0(low);
            this->createPleaveFactorTable(pleaveFactorTable);
            this->createPleaveaTable(pleaveaTable, pleaveFactorTable);

            
            const Real low_value(GSL_FN_EVAL(&F, low));
            
            if (low_value <= 0.0)
            {
                break;
            }
            
            // FIXME: 
            if (fabs(low) <= minT || 
                fabs(low_value - low_value_prev) < TOLERANCE) 
            {
                log_.info("couldn't adjust low. Fa(%.16g) = %.16g; r0 = %.16g, %s",
                           low, GSL_FN_EVAL(&F, low), r0, dump().c_str());
                log_.info("returning %.16g", minT);
                return minT;
            }
            low_value_prev = low_value;

            log_.info("drawTime2: adjusting low: %.16g, Fa = %.16g", low, low_value);
            low *= .1;
        }
    }

    const Real t(findRoot(F, solver, low, high, 0.,
                            this->TOLERANCE, "drawTime2: a"));

    return t;
}


Real 
GreensFunction3DRadAbs::drawPleaves(gsl_function const& F,
                                            gsl_root_fsolver* solver,
                                            Real t_guess,
                                            RealVector& pleaveFactorTable,
                                            RealVector& pleavesTable) const
{
    const Real minT(1e-12);

    Real low(t_guess);
    Real high(t_guess);

    // adjust high and low to make sure that f(low) and f(high) straddle.
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

            if (fabs(high) >= 1e10)
            {
                throw std::runtime_error(
                    (boost::format(
                        "couldn't adjust high. Fs(%.16g) = %.16g; r0 = %.16g, %s") %
                        high % GSL_FN_EVAL(&F, high) % r0 % dump()).str());
            }

            log_.info("drawTime2: adjusting high: %.16g Fs = %.16g", 
                      high, high_value);
            high *= 10;
        }
    }
    else
    {
        Real low_value_prev(value);
        low *= .1;

        for (;;)
        {
            this->updateAlphaTable0(low);
            this->createPleaveFactorTable(pleaveFactorTable);
            this->createPleavesTable(pleavesTable, pleaveFactorTable);
            
            const Real low_value(GSL_FN_EVAL(&F, low));
            
            if (low_value <= 0.0)
            {
                break;
            }
            
            // FIXME: 
            if (fabs(low) <= minT)//|| 
//                fabs(low_value - low_value_prev) < TOLERANCE) 
            {
                log_.info("couldn't adjust low.  returning minT (=%.16g);"
                          "Fs(%.16g) = %.16g; r0 = %.16g, %s",
                          minT, low, GSL_FN_EVAL(&F, low), r0, dump().c_str());
                return minT;
            }
            low_value_prev = low_value;

            log_.info("drawTime2: adjusting low: %.16g, Fs = %.16g", low, low_value);
            low *= .1;
        }
    }

    const Real t(findRoot(F, solver, low, high, 0., this->TOLERANCE,
                            "drawTime2: s"));

    return t;
}




Real GreensFunction3DRadAbs::drawR(Real rnd, Real t) const
{
    const Real D(this->getD());
    const Real sigma(this->getSigma());
    const Real a(this->geta());

    if (!(rnd < 1.0 && rnd >= 0.0))
    {
        throw std::invalid_argument((boost::format("rnd < 1.0 && rnd >= 0.0 : rnd=%.16g") % rnd).str());
    }

    if (!(r0 >= sigma && r0 < a))
    {
        throw std::invalid_argument((boost::format("r0 >= sigma && r0 < a : r0=%.16g, sigma=%.16g, a=%.16g") % r0 % sigma % a).str());
    }


    if (t == 0.0)
    {
        return r0;
    }

    const Real psurv(p_survival(t));

//    RealVector num_r0Table;
//    createNum_r0Table(num_r0Table, r0);

    p_int_r_params params = { this, t, /*num_r0Table,*/ rnd * psurv };

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
    if (GSL_FN_EVAL(&F, r0) < 0.0)
    {
        // low = r0
        unsigned int H(3);

        for (;;)
        {
            high = r0 + H * sqrt6Dt;
            if (high > a)
            {
                if (GSL_FN_EVAL(&F, a) < 0.0)
                {
                    log_.info("drawR: p_int_r(a) < 0.0. returning a");
                    return a;
                }

                high = a;
                break;
            }

            const Real value(GSL_FN_EVAL(&F, high));
            if (value > 0.0)
            {
                break;
            }

            ++H;
        }

    }
    else
    {
        // high = r0
        unsigned int H(3);

        for (;;)
        {
            low = r0 - H * sqrt6Dt;
            if (low < sigma)
            {
                if (GSL_FN_EVAL(&F, sigma) > 0.0)
                {
                    log_.info("drawR: p_int_r(sigma) > 0.0. returning sigma");
                    return sigma;
                }

                low = sigma;
                break;
            }

            const Real value(GSL_FN_EVAL(&F, low));
            if (value < 0.0)
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




Real GreensFunction3DRadAbs::p_n_alpha(unsigned int i, unsigned int n,
                                               Real r, Real t) const
{
    const Real sigma(this->getSigma());
    const Real h(this->geth());
    const Real a(this->geta());

    const Real mDt(- this->getD() * t);

    const Real alpha(this->getAlpha(n, i));
    const Real alphasq(alpha * alpha);

    const Real aAlpha(a * alpha);
    const Real sigmaAlpha(sigma * alpha);
    const Real hSigma(geth() * getSigma());
    const Real realn(static_cast<Real>(n));
    const Real hSigma_m_n(hSigma - realn);

    const Real term1(alphasq * alphasq * exp(mDt * alphasq));


    const SphericalBesselGenerator& s(SphericalBesselGenerator::instance());

    const Real js1(s.j(n,   sigmaAlpha));
    const Real js2(s.j(n+1, sigmaAlpha));
    const Real ja( s.j(n,   aAlpha));
    const Real ya( s.y(n,   aAlpha));
    const Real jr( s.j(n,   r * alpha));
    const Real yr( s.y(n,   r * alpha));
    const Real jr0(s.j(n,   r0 * alpha));
    const Real yr0(s.y(n,   r0 * alpha));

    const Real J(hSigma_m_n * js1 + sigmaAlpha * js2);
    const Real Jsq(J * J);

    const Real JY1(ja * yr - ya * jr);
    const Real JY2(ja * yr0 - ya * jr0);

    const Real num(Jsq * JY1 * JY2);

    const Real den1(a * (realn + realn * realn - 
                           sigma * (h + h * h * sigma + sigma * alphasq))
                     * ja * ja);

    const Real den2(sigma * Jsq);

    const Real den(den1 + den2);

    const Real result(term1 * num / den);

    return result;
}



Real 
GreensFunction3DRadAbs::p_n(Integer n, Real r, Real t, Real max_alpha) const
{
    const unsigned int min_i(2);

    Real p(0.0);
    
    Integer i(0);
    for (;;)
    {
        const Real alpha(getAlpha(n,i));

        const Real p_i(p_n_alpha(i, n, r, t));
        p += p_i;

        if(alpha >= max_alpha && i >= min_i)
        {
            break;
        }

        if(i == MAX_ALPHA_SEQ)
        {
            break;
        }

        ++i;
    }

    return p;
}

void
GreensFunction3DRadAbs::makep_nTable(RealVector& p_nTable,
                                             Real r, Real t) const
{
    const Real sigma(this->getSigma());
    const Real a(this->geta());

    p_nTable.clear();

    const Real factor(a * sigma / (M_PI * 2));

    const Real Dt(this->getD() * t);
    const Real alpha00(this->getAlpha(0, 0));

    const Real max_alpha(sqrt(alpha00 * alpha00 - 
                              log(THETA_TOLERANCE * 1e-1) / Dt));


    const Real p_0(this->p_n(0, r, t, max_alpha) * factor);

    p_nTable.push_back(p_0);

    if(p_0 == 0)
    {
        return;
    }

    const Real threshold(fabs(THETA_TOLERANCE * p_0));

    Real p_n_prev_abs(fabs(p_0));
    unsigned int n(1);
    for (;;)
    {
        if(getAlpha(n, 0) >= max_alpha)
        {
            break;
        }

        Real p_n(this->p_n(n, r, t, max_alpha) * factor);
        
        p_nTable.push_back(p_n);
        const Real p_n_abs(fabs(p_n));
        // truncate when converged enough.
        if(p_n_abs < threshold &&
           p_n_prev_abs < threshold &&
           p_n_abs <= p_n_prev_abs)
        {
            break;
        }
        
        if (n >= this->MAX_ORDER)
        {
            break;
        }
        
        ++n;
        p_n_prev_abs = p_n_abs;
    }
}


Real 
GreensFunction3DRadAbs::dp_n_alpha_at_a(unsigned int i, unsigned int n,
                                                Real t) const
{
    const Real sigma(this->getSigma());
    const Real h(this->geth());
    const Real a(this->geta());

    const Real mDt(- this->getD() * t);
    const Real alpha(this->getAlpha(n, i));

    const Real alphasq(alpha * alpha);

    const Real aAlpha(a * alpha);
    const Real sigmaAlpha(sigma * alpha);
    const Real hSigma(geth() * getSigma());
    const Real realn(static_cast<Real>(n));
    const Real hSigma_m_n(hSigma - realn);

    const Real term1(alphasq * alpha * exp(mDt * alphasq));

    const SphericalBesselGenerator& s(SphericalBesselGenerator::instance());

    const Real js1(s.j(n,   sigmaAlpha));
    const Real js2(s.j(n+1, sigmaAlpha));
    const Real ja( s.j(n,   aAlpha));
    const Real ya( s.y(n,   aAlpha));
    const Real jr0(s.j(n,   r0 * alpha));
    const Real yr0(s.y(n,   r0 * alpha));

    const Real J(hSigma_m_n * js1 + sigmaAlpha * js2);
    const Real Jsq(J * J);

    const Real JY(- jr0 * ya + ja * yr0);

    const Real num(Jsq * JY);

    const Real den1(a * (realn + realn * realn - 
                           sigma * (h + h * h * sigma + sigma * alphasq))
                     * ja * ja);

    const Real den2(sigma * Jsq);

    const Real den(den1 + den2);

    const Real result(term1 * num / den);

    return result;
}

Real 
GreensFunction3DRadAbs::dp_n_at_a(Integer n, Real t,
                                          Real max_alpha) const
{
    const unsigned int min_i(2);

    Real p(0.0);
    
    Integer i(0);
    for (;;)
    {
        const Real alpha(getAlpha(n,i));

        const Real p_i(dp_n_alpha_at_a(i, n, t));

        p += p_i;

        if(alpha >= max_alpha && i >= min_i)
        {
            break;
        }

        if(i == MAX_ALPHA_SEQ)
        {
            break;
        }

        ++i;
    }

    return p;
}


void
GreensFunction3DRadAbs::makedp_n_at_aTable(RealVector& p_nTable,
                                                   Real t) const
{
    const Real sigma(this->getSigma());
    const Real a(this->geta());

    p_nTable.clear();

    const Real factor(this->getD() * sigma / (a * M_PI * 2));

    const Real Dt(this->getD() * t);
    const Real alpha00(this->getAlpha(0, 0));

    const Real max_alpha(sqrt(Dt * alpha00 * alpha00 - 
                              log(THETA_TOLERANCE * 1e-1) / Dt));


    const Real p_0(this->dp_n_at_a(0, t, max_alpha) * factor);

    p_nTable.push_back(p_0);

    if(p_0 == 0)
    {
        return;
    }

    const Real threshold(fabs(THETA_TOLERANCE * p_0));

    Real p_n_prev_abs(fabs(p_0));
    unsigned int n(1);
    for (;;)
    {
        if(getAlpha(n, 0) >= max_alpha)
        {
            break;
        }

        Real p_n(this->dp_n_at_a(n, t, max_alpha) * factor);

        p_nTable.push_back(p_n);
        
        const Real p_n_abs(fabs(p_n));
        // truncate when converged enough.
        if(p_n_abs < threshold &&
           p_n_prev_abs < threshold &&
           p_n_abs <= p_n_prev_abs)
        {
            break;
        }
        
        if (n >= this->MAX_ORDER)
        {
            break;
        }
        
        ++n;
        p_n_prev_abs = p_n_abs;
    }

}

Real 
GreensFunction3DRadAbs::p_theta(Real theta, Real r, Real t) const 
{
    {
        const Real sigma(this->getSigma());
        const Real a(this->geta());
        
        if (!(theta >= 0.0 && theta <= M_PI))
        {
            throw std::invalid_argument((boost::format("theta >= 0.0 && theta <= M_PI : theta=%.16g, M_PI=%.16g") % theta % M_PI).str());
        }

        // r \in (sigma, a);  not defined at r == sigma and r == a.
        if (!(r >= sigma && r < a))
        {
            throw std::invalid_argument((boost::format("r >= sigma && r < a : r=%.16g, sigma=%.16g, a=%.16g") % r % sigma % a).str());
        }

        if (!(r0 >= sigma && r0 < a))
        {
            throw std::invalid_argument((boost::format("r0 >= sigma && r0 < a : r0=%.16g, sigma=%.16g, a=%.16g") % r0 % sigma % a).str());
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

Real GreensFunction3DRadAbs::dp_theta(Real theta, Real r, Real t) const 
{
    {
        const Real sigma(this->getSigma());
        const Real a(this->geta());
        
        if (!(theta >= 0.0 && theta <= M_PI))
        {
            throw std::invalid_argument((boost::format("theta >= 0.0 && theta <= M_PI : theta=%.16g, M_PI=%.16g") % theta % M_PI).str());
        }


        // r \in [ sigma, a ]  ;  unlike p_theta,
        // defined at r == sigma and r == a.
        if (!(r >= sigma && r <= a))
        {
            throw std::invalid_argument((boost::format("r >= sigma && r <= a : r=%.16g, sigma=%.16g, a=%.16g") % r % sigma % a).str());
        }

        if (!(r0 >= sigma && r0 < a))
        {
            throw std::invalid_argument((boost::format("r0 >= sigma && r0 < a : r0=%.16g, sigma=%.16g, a=%.16g") % r0 % sigma % a).str());
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

    makedp_n_at_aTable(p_nTable, t);

    const Real p(p_theta_table(theta, r, t, p_nTable));

    return p;
}


static Real
p_theta_n(unsigned int n,
          GreensFunction3DRadAbs::RealVector const& p_nTable,
          GreensFunction3DRadAbs::RealVector const& lgndTable)
{
    return p_nTable[n] * lgndTable[n] * (2 * n + 1);
}

Real
GreensFunction3DRadAbs::p_theta_table(Real theta, Real r,
                                              Real t,
                                              RealVector const& p_nTable) const
{
    const unsigned int tableSize(p_nTable.size());

    Real sin_theta;
    Real cos_theta;
    sincos(theta, &sin_theta, &cos_theta);

    RealVector lgndTable(tableSize);
    gsl_sf_legendre_Pl_array(tableSize-1, cos_theta, &lgndTable[0]);

    return funcSum_all(
            boost::bind(&p_theta_n, _1, p_nTable, lgndTable),
            tableSize) * sin_theta;
}

void
GreensFunction3DRadAbs::
make_p_thetaTable(RealVector& pTable,
                  Real r, 
                  Real t,
                  unsigned int n,
                  RealVector const& p_nTable) const
{
    const Real thetaStep(M_PI / n);

    pTable.push_back(0.0);

    Real p_prev(0.0);
    unsigned int i(1);
    for (;;)
    {
        const Real theta(thetaStep * i);

        Real p(this->p_theta_table(theta, r, t, p_nTable));

        if (p < 0.0)
        {
            log_.info("drawTheta: p<0 %.16g", p);
            p = 0.0;
        }

        const Real value((p_prev + p) * 0.5);
        pTable.push_back(*(pTable.end() - 1) + value);

        if (/* value < pTable[i] * std::numeric_limits<Real>::epsilon() || */
            i >= n - 1)
        {
            break;   // pTable is valid in [0,i].
        }

        p_prev = p;
        ++i;
    }

}


Real 
GreensFunction3DRadAbs::ip_theta(Real theta, Real r, Real t) const
{
    {
        const Real sigma(this->getSigma());
        const Real a(this->geta());
        
        if (!(theta >= 0.0 && theta <= M_PI))
        {
            throw std::invalid_argument((boost::format("theta >= 0.0 && theta <= M_PI : theta=%.16g, M_PI=%.16g") % theta % M_PI).str());
        }

        // r \in (sigma, a)
        if (!(r >= sigma && r < a))
        {
            throw std::invalid_argument((boost::format("r >= sigma && r < a : r=%.16g, sigma=%.16g, a=%.16g") % r % sigma % a).str());
        }

        if (!(r0 >= sigma && r0 < a))
        {
            throw std::invalid_argument((boost::format("r0 >= sigma && r0 < a : r0=%.16g, sigma=%.16g, a=%.16g") % r0 % sigma % a).str());
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


Real 
GreensFunction3DRadAbs::idp_theta(Real theta, Real r, Real t) const
{
    {
        const Real sigma(this->getSigma());
        const Real a(this->geta());
        
        if (!(theta >= 0.0 && theta <= M_PI))
        {
            throw std::invalid_argument((boost::format("theta >= 0.0 && theta <= M_PI : theta=%.16g, M_PI=%.16g") % theta % M_PI).str());
        }

        // r \in [ sigma, a ]
        if (!(r >= sigma && r <= a))
        {
            throw std::invalid_argument((boost::format("r >= sigma && r <= a : r=%.16g, sigma=%.16g, a=%.16g") % r % sigma % a).str());
        }

        if (!(r0 >= sigma && r0 < a))
        {
            throw std::invalid_argument((boost::format("r0 >= sigma && r0 < a : r0=%.16g, sigma=%.16g, a=%.16g") % r0 % sigma % a).str());
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

    makedp_n_at_aTable(p_nTable, t);

    const Real p(ip_theta_table(theta, r, t, p_nTable));

    return p;
}

static Real
ip_theta_n(unsigned int n,
           GreensFunction3DRadAbs::RealVector const& p_nTable,
           GreensFunction3DRadAbs::RealVector const& lgndTable1)
{
    // lgndTable1 is offset by 1; lgndTable1[0] is for n=-1.

    const Real lgnd_n_m1(lgndTable1[n]);   // n-1
    const Real lgnd_n_p1(lgndTable1[n+2]); // n+1
    
    // the term (1 + 2 n) is canceled out.
    return p_nTable[n] * (lgnd_n_m1 - lgnd_n_p1);
}


Real 
GreensFunction3DRadAbs::ip_theta_table(Real theta, Real r,
                                               Real t, RealVector const& p_nTable) const
{
    const unsigned int tableSize(p_nTable.size());

    const Real cos_theta(cos(theta));

    // LgndTable is offset by 1 to incorporate the n=-1 case.
    // For ex: LgndTable[0] is for n=-1, lgndTable[1] is for n=0 ...

    RealVector lgndTable(tableSize + 2);
    lgndTable[0] = 1.0;  // n = -1
    gsl_sf_legendre_Pl_array(tableSize, cos_theta, &lgndTable[1]);

    return funcSum_all(
        boost::bind(&ip_theta_n, _1, p_nTable, lgndTable),
        tableSize);
}

struct GreensFunction3DRadAbs::ip_theta_params
{ 
    GreensFunction3DRadAbs const* const gf;
    const Real r;
    const Real t;
    RealVector const& p_nTable;
    const Real value;
};

Real GreensFunction3DRadAbs::ip_theta_F(Real theta, ip_theta_params const* params)
{
    const GreensFunction3DRadAbs* const gf(params->gf); 
    const Real r(params->r);
    const Real t(params->t);
    RealVector const& p_nTable(params->p_nTable);
    const Real value(params->value);

    return gf->ip_theta_table(theta, r, t, p_nTable) - value;
}

Real 
GreensFunction3DRadAbs::drawTheta(Real rnd, Real r, Real t) const
{
    Real theta;

    const Real sigma(this->getSigma());
    const Real a(this->geta());

    // input parameter range checks.
    if (!(rnd < 1.0 && rnd >= 0.0))
    {
        throw std::invalid_argument((boost::format("rnd < 1.0 && rnd >= 0.0 : rnd=%.16g") % rnd).str());
    }

    if (!(r0 >= sigma && r0 < a))
    {
        throw std::invalid_argument((boost::format("r0 >= sigma && r0 < a : r0=%.16g, sigma=%.16g, a=%.16g") % r0 % sigma % a).str());
    }

    if (!(r >= sigma))
    {
        throw std::invalid_argument((boost::format("r >= sigma : r=%.16g, sigma=%.16g") % r % sigma).str());
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

    const Real high(M_PI);

    RealVector p_nTable;

    if (r >= geta())
    {
        //puts("dp");
        makedp_n_at_aTable(p_nTable, t);
    }
    else
    {
        makep_nTable(p_nTable, r, t);
    }

    const Real ip_theta_pi(ip_theta_table(high, r, t, p_nTable));

    ip_theta_params params = { this, r, t, p_nTable, rnd * ip_theta_pi };

    gsl_function F = 
        { reinterpret_cast<typeof(F.function)>(&ip_theta_F), &params };

    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));
    gsl_root_fsolver_set(solver, &F, 0.0, high);

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


//
// debug
//

std::string GreensFunction3DRadAbs::dump() const
{
    std::ostringstream ss;
    ss << "D = " << this->getD() << 
        ", r0 = " << this->getr0() <<
        ", sigma = " << this->getSigma() <<
        ", a = " << this->geta() <<
        ", kf = " << this->getkf() <<
        ", h = " << this->geth() << std::endl;
    return ss.str();
}

Logger& GreensFunction3DRadAbs::log_(
        Logger::get_logger("GreensFunction3DRadAbs"));

#if !defined( __FIRSTPASSAGEGREENSFUNCTION1D_HPP )
#define __FIRSTPASSAGEGREENSFUNCTION1D_HPP

#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_roots.h>

#include <math.h>

#include "findRoot.hpp"
#include "Defs.hpp"
#include "OldDefs.hpp"			// TODO: this must be removed at some point!

#include "GreensFunction.hpp"
#include "PairGreensFunction.hpp"	// needed to declare EventType


class GreensFunction1DAbsAbs: public GreensFunction
{
private:
    // This is a typical length scale of the system, may not be true!
    static const Real L_TYPICAL = 1E-8;
    // The typical timescale of the system, may also not be true!!
    static const Real T_TYPICAL = 1E-6;
    // measure of 'sameness' when comparing floating points numbers
    static const Real EPSILON = 1E-12;
    //E3; Is 1E3 a good measure for the probability density?!
    static const Real PDENS_TYPICAL = 1;
    // The maximum number of terms in the sum
    static const int MAX_TERMS = 500;
    // The minimum
    static const int MIN_TERMS = 20;

public:
    enum EventKind
    {
        IV_ESCAPE,
        IV_REACTION
    };

public:
    GreensFunction1DAbsAbs(Real D, Real r0, Real sigma, Real a)
	: GreensFunction(D), v(0.0), sigma(sigma), a(a), r0(r0), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
    {
	;   // do nothing
    }

    // The constructor is overloaded and can be called with or without drift v
    GreensFunction1DAbsAbs(Real D, Real v, Real r0, Real sigma, Real a) // copy constructor including drift variable v
	: GreensFunction(D), v(v), sigma(sigma), a(a), r0(r0), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
    {
	;   // do nothing
    }

    ~GreensFunction1DAbsAbs()
    { 
	;   // empty
    }

    // This also sets the scale
    void seta(Real a)
    {
	Real L( a - this->sigma );
	
	THROW_UNLESS( std::invalid_argument, L >= 0.0 && (this->r0 - sigma) <= L);

	// Use a typical domain size to determine if we are here
	// defining a domain of size 0.
	if ( L <= EPSILON * l_scale )
	{
	    // just some random value to show that the domain is zero
	    this->a = -INT_MAX;
	}
	else
	{   
	    // set the typical time scale (msd = sqrt(2*d*D*t) )
	    // this is needed by drawTime_f, do not get rid of it!
	    this->t_scale = (L*L)/this->getD();
	    // set a
	    this->a = a;
	}
    }

    Real getsigma() const
    {
	return this->sigma;
    }
 
    Real geta() const
    {
	return this->a;
    }

    Real getv() const
    {
	return this->v;
    }

    void setr0(Real r0)
    {
	if ( this->a - this->sigma < 0.0 )
	{
	    // if the domain had zero size    
	    THROW_UNLESS( std::invalid_argument,
	                  0.0 <= (r0-sigma) && (r0-sigma) <= EPSILON * l_scale );
	    this->r0 = 0.0;
	}
	else
	{
	    // The normal case
	    THROW_UNLESS( std::invalid_argument,
	                  0.0 <= (r0-sigma) && r0 <= this->a);
	    this->r0 = r0;
	}
    }

    Real getr0() const
    {
	return this->r0;
    }

    // Draws the first passage time from the propensity function
    Real drawTime (Real rnd) const;

    // Draws the position of the particle at a given time, assuming that 
    // the particle is still in the
    // domain
    Real drawR (Real rnd, Real t) const;

    // Calculates the amount of flux leaving the left boundary at time t
    Real leaves(Real t) const;

    // Calculates the amount of flux leaving the right boundary at time t
    Real leavea(Real t) const;

    // Determines based on the flux ratios if the particle left the left 
    // or right boundary
    EventKind drawEventType( Real rnd, Real t ) const;

    // Calculates the probability of finding the particle inside the 
    // domain at time t so, the survival probability
    Real p_survival (Real t) const;

    // Calculates the probability density of finding the particle at 
    // location z at timepoint t, given that the particle is still in the 
    // domain.
    Real calcpcum (Real r, Real t) const;
    
    // Calculates the probability density of finding the particle at 
    // location r at time t.
    Real prob_r (Real r, Real t) const;

    std::string dump() const;

private:
    struct drawT_params
    {
	// use 10 terms in the summation for now
	double exponent[MAX_TERMS];
	double Xn[MAX_TERMS];
	double prefactor;
	int    terms;
	Real tscale;
	// random number
	double rnd;
    };

    static double drawT_f (double t, void *p);

    struct drawR_params
    {
	double S_Cn_An[MAX_TERMS];
	double n_L[MAX_TERMS];
	// variables H: for additional terms appearing as multiplicative factors etc.
	double H[5];
	int terms;
	// random number
	double rnd;
    };

    static double drawR_f (double z, void *p);

private:
    // The diffusion constant and drift velocity
    Real v;
    // These are the dimensions of our domain; L is calculated as a-sigma
    Real sigma;
    Real a;
    Real r0;
    // This is the 'length scale' of your system (1e-14 or 1e6)
    // Although rescaling is discontinued, we use it to check whether a is well-chosen
    Real l_scale;
    // This is the time scale of the system, used by drawTime_f
    Real t_scale;
};
#endif // __FIRSTPASSAGEGREENSFUNCTION1D_HPP

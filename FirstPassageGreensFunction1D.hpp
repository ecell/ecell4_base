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
#include "PairGreensFunction.hpp"	// needed to declare EventType
					// This is new since being ported to the new version of the code
					// -> TODO: Improve organisation of the included files
					//
					// The ESCAPE and REACTION events have been changed to IV_ESCAPE / IV_REACTION

class FirstPassageGreensFunction1D
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
    static const int MAX_TERMEN = 500;
    // The minimum
    static const int MIN_TERMEN = 20;

public:
    FirstPassageGreensFunction1D(Real D, Real r0, Real a)	// standard version as used up to now
	: D(D), v(0.0), a(a), r0(r0), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
    {
	;   // do nothing
    }

    // The constructor is overloaded and can be called with or without drift v
    FirstPassageGreensFunction1D(Real D, Real v, Real r0, Real a) // copy constructor including drift variable v
	: D(D), v(v), a(a), r0(r0), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
    {
	;   // do nothing
    }

    ~FirstPassageGreensFunction1D()
    { 
	;   // empty
    }

    // This also sets the scale
    void seta(Real a)
    {
	THROW_UNLESS( std::invalid_argument, a >= 0.0 && r0 <= a);

	// Use a typical domain size to determine if we are here
	// defining a domain of size 0.
	if ( a <= EPSILON * l_scale )
	{
	    // just some random value to show that the domain is 
	    // zero
	    this->a = -1.0;
	    // don't touch the scales
	    //this->l_scale = 1.0;
	}
	else
	{
	    // set the scale to the given one
	    //this->l_scale = a;		// renormalized version, discontinued
	    // set the typical time scale (msd = sqrt(2*d*D*t) )
	    this->t_scale = (a*a)/D;		// this is needed by drawTime_f, do not get rid of it!
	    // set a
	    this->a = a;
	    // this->a = 1.0;			// renormalized version, discontinued
	}
    }

    Real geta() const
    {
	return this->a;
    }

    Real getD() const
    {
	return this->D;
	// return this->D/(l_scale*l_scale);	// renormalized version, discontinued
    }

    Real getv() const
    {
	return this->v;
	// return this->v/l_scale;		// renormalized version, discontinued
    }

    void setr0(Real r0)
    {
	if ( this->a < 0.0 )
	{
	    // if the domain had zero size    
	    THROW_UNLESS( std::invalid_argument,
	                  0.0 <= r0 && r0 <= EPSILON * l_scale );
	    this->r0 = 0.0;
	}
	else
	{
	    // The normal case
	    THROW_UNLESS( std::invalid_argument,
	                  0.0 <= r0 && r0 <= this->a);
	    this->r0 = r0;
	}
    }

    Real getr0() const
    {
	return this->r0;
	// return this->r0/l_scale;		// renormalized version, discontinued
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
    EventType drawEventType( Real rnd, Real t ) const;

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

private:
    struct drawT_params
    {
	// use 10 terms in the summation for now
	double exponent[MAX_TERMEN];
	double Xn[MAX_TERMEN];
	double prefactor;
	int    terms;
	Real tscale;
	// the random number associated with the time
	double rnd;
    };

    static double drawT_f (double t, void *p);

    struct drawR_params
    {
	double S_Cn_An[MAX_TERMEN];
	double n_l[MAX_TERMEN];
	// constants needed for the sum computation in case of nonzero drift
	double v2D;
	int terms;
	// the random number associated with the time
	double rnd;
    };

    static double drawR_f (double z, void *p);

private:
    // The diffusion constant and drift velocity
    Real D;
    Real v;
    // The length of your domain (also the l_scale, see below)
    Real a;
    Real r0;
    // This is the 'length scale' of your system (1e-14 or 1e6)
    // Although rescaling is discontinued, we use it to check whether a is well-chosen
    Real l_scale;
    // This is the time scale of the system, used by drawTime_f
    Real t_scale;
};
#endif // __FIRSTPASSAGEGREENSFUNCTION1D_HPP

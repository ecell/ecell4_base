#if !defined( __FIRSTPASSAGEGREENSFUNCTION1DRAD_HPP )
#define __FIRSTPASSAGEGREENSFUNCTION1DRAD_HPP

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

class FirstPassageGreensFunction1DRad: public GreensFunction
{
private:
    // This is a typical length scale of the system, may not be true!
    static const Real L_TYPICAL = 1E-8;
    // The typical timescale of the system, may also not be true!!
    static const Real T_TYPICAL = 1E-6;
    // measure of 'sameness' when comparing floating points numbers
    static const Real EPSILON = 1E-10;
    // Is 1E3 a good measure for the probability density?!
    static const Real PDENS_TYPICAL = 1;
    // The maximum number of terms used in calculating the sum
    static const int MAX_TERMEN = 500;
    // The minimum number of terms
    static const int MIN_TERMEN = 20;


public:
    FirstPassageGreensFunction1DRad(Real D, Real k, Real r0, Real sigma, Real a)
	: GreensFunction(D), v(0.0), k(k), r0(r0), sigma(sigma), a(a), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
    {
	// do nothing
    }

    // The constructor is overloaded and can be called with or without drift v
    // copy constructor including drift variable v
    FirstPassageGreensFunction1DRad(Real D, Real k, Real v, Real r0, Real sigma, Real a)
	: GreensFunction(D), v(v), k(k), r0(r0), sigma(sigma), a(a), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
    {
	// do nothing
    }

    ~FirstPassageGreensFunction1DRad()
    {
	;   // empty
    }

    // This also sets the scale
    void seta(Real a)
    {
	THROW_UNLESS( std::invalid_argument, (a-this->sigma) >= 0.0 && this->r0 <= a);

	// Use a typical domain size to determine if we are here 
	// defining a domain of size 0.
	if ( (a-this->sigma) < EPSILON * this->l_scale)
	{
	    // just some random value to show that the domain is 
	    // zero
	    this->a = -1.0;
	    //this->l_scale = 1.0;
	}
	else
	{
	    // set the l_scale to the given one
	    this->l_scale = a-sigma;
	    // set the typical time scale (msd = sqrt(2*d*D*t) )
	    this->t_scale = (l_scale*l_scale)/this->getD();
	    // this->a = a/l_scale		// renormalized version, discontinued
	    this->a = a;
	}
    }

    Real geta() const
    {
	return this->a;
    }
    
    Real getsigma() const
    {
	return this->sigma;
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
	return r0;
	//return r0/l_scale;				// renormalized version, discontinued
    }

    Real getk() const
    {
	return this->k;
	// don't forget to scale the k as well!
	//return this->k/l_scale;			// renormalized version, discontinued
    }

    Real getv() const
    {
	return this->v;
	//return this->v/l_scale;			// renormalized version, discontinued
    }

    // Calculates the probability density of finding the particle at 
    // location z at timepoint t, given that the particle is still in the 
    // domain.
    Real calcpcum (Real r, Real t) const;

    // Determine which event has occured, an escape or a reaction. Based 
    // on the fluxes through the boundaries at the given time. Beware: if 
    // t is not a first passage time you still get an answer!
    EventType drawEventType( Real rnd, Real t ) const;

    // Draws the first passage time from the propensity function
    Real drawTime (Real rnd) const;

    // Draws the position of the particle at a given time, assuming that 
    // the particle is still in the domain
    Real drawR (Real rnd, Real t) const;


// These methods are both public and private, they are used by public methods 
// but can also be called from the 'outside'. This is mainly because of 
// debugging purposes.


    // Calculates the probability of finding the particle inside the 
    // domain at time t -> the survival probability
    Real p_survival (Real t) const;

    // Calculates the total probability flux leaving the domain at time t
    Real flux_tot (Real t) const;

    // Calculates the probability flux leaving the domain through the 
    // radiative boundary at time t
    Real flux_rad (Real t) const;

    // Calculates the flux leaving the domain through the radiative 
    // boundary as a fraction of the total flux. This is the probability 
    // that the particle left the domain through the radiative
    // boundary instead of the absorbing boundary.
    Real fluxRatioRadTot (Real t) const;

    // Calculates the probability density of finding the particle at 
    // location r at time t.
    Real prob_r (Real r, Real t) const;
    
// End of public/private mix methods

//private:
    // Calculates the roots of tan(a*x)=-xk/h
    Real root_n(int n) const;
    
private:

    Real An (Real a_n) const;

    Real Bn (Real a_n) const;

    Real Cn (Real a_n, Real t) const;

    struct tan_f_params
    {
	Real a;
	Real h;
    };

    static double tan_f (double x, void *p);
    // this is the appropriate definition of the function in gsl

    struct drawT_params
    {
	double exponent[MAX_TERMEN];
	double Xn[MAX_TERMEN];
	double prefactor;
	int    terms;
	// the timescale used for convergence
	Real   tscale;
	// the random number associated with the time
	double rnd;
    };

    static double drawT_f (double t, void *p);

    struct drawR_params
    {
	double root_n[MAX_TERMEN];
	double S_Cn_root_n[MAX_TERMEN];
	// variables H: for additional terms appearing as multiplicative factors etc.
	double H[5];
	int terms;
	// the random number associated with the time
	double rnd;
    };

    static double drawR_f (double z, void *p);

    // The diffusion constant and drift velocity
    Real v;
    // The reaction constant
    Real k;
    Real r0;
    // The left and right boundary of the domain (also the l_scale, see below)
    Real sigma;
    Real a;
    // This is the 'scale' of the system (1e-14 or 1e6).
    Real l_scale;
    // This is the time scale of the system.
    Real t_scale;
};
#endif // __FIRSTPASSAGEGREENSFUNCTION1DRAD_HPP

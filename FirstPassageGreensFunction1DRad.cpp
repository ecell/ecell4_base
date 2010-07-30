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

#include "FirstPassageGreensFunction1DRad.hpp"

// this is the appropriate definition of the function in gsl
double FirstPassageGreensFunction1DRad::tan_f (double x, void *p)
{
    // casts the void to the struct pointer
    struct tan_f_params *params = (struct tan_f_params *)p;
    const Real a = (params->a);
    const Real h = (params->h);
    const Real h_a (h*a);
    if ( h_a < 1 )
    {
	// h = k/D, h=h1/k1
	return 1/tan(x) + (h_a)/x;
    }
    else
    {
	// h = k/D, h=h1/k1
	return tan(x) + x/(h_a);
    }
}

// Calculates the roots of tan(x*a)=-xk/h
Real FirstPassageGreensFunction1DRad::a_n(int n) const
{
    // a=length of domain
    const Real a(this->geta());
    // h=(k+v/2)/D
    const Real h((this->getk()+this->getv()/2.0)/this->getD());	// the drift v also comes into this constant
    Real upper, lower;

    if ( h*a < 1 )
    {
	// 1E-10 to make sure that he doesn't include the transition 
	lower = (n-1)*M_PI + 1E-10;
	// (asymptotic) from infinity to -infinity
	upper =  n   *M_PI - 1E-10;
    }
    else
    {
	lower = (n-1)*M_PI + M_PI_2 + 1E-10;
	upper = n    *M_PI + M_PI_2 - 1E-10;
    }

    gsl_function F;
    struct tan_f_params params = { a, h };
     
    F.function = &FirstPassageGreensFunction1DRad::tan_f;
    F.params = &params;

    // define a new solver type brent
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );

    // make a new solver instance
    // TODO: incl typecast?
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    // get the root = run the rootsolver
    const Real root( findRoot( F, solver, lower, upper, 1.0*EPSILON, EPSILON,
                            "FirstPassageGreensFunction1DRad::root_tan" ) );
    gsl_root_fsolver_free( solver );
    
    return root/a;
    // This rescaling is important, because the function tan_f is used to solve for
    // tan(x)+x/h/a=0, whereas we actually need tan(x*a)+x/h=0, So if x solves the 
    // subsidiary equation, x/a solves the original one.
}

// This is the non-exponential factor in the Green's function sum, not
// including the factor containing the explicit r-dependency (The latter
// is given by the Bn's, see below).
//
// r0 is here still in the interval from 0 to a (and supposed to be the
// starting point of the particle at t0).
//
// The root a_n also must be the specific one for that interval, thus
// the one rescaled by a (see comments in function a_n(n) ).
//
// The factor calculated here is identical for the cases w. or w/o drift,
// only h changes.
Real FirstPassageGreensFunction1DRad::An (Real a_n) const
{
    const Real h((this->getk()+this->getv()/2.0)/this->getD());
    const Real a(this->geta());
    const Real r0(this->getr0());
    const Real anr0 = a_n*r0;

    return (a_n*cos(anr0) + h*sin(anr0)) / (h + (a_n*a_n + h*h)*a);
}

// r is here is in the interval [0,a]
Real FirstPassageGreensFunction1DRad::Bn (Real a_n) const
{
    const Real h((this->getk()+this->getv()/2.0)/this->getD());
    const Real k(this->getk());
    const Real a(this->geta());
    const Real D(this->getD());
    const Real v(this->getv());
    
    const Real ana(a_n*a);
    const Real an2(a_n*a_n);
    const Real h2(h*h);
    const Real v2D(v/2.0/D);

    if(v==0.0)	return (h2 - (an2 + h2)*cos(ana)) / (h*a_n);
    else	return (h*k/D - exp(v2D*a)*(an2+h2)*cos(ana) ) / (h/a_n*(an2+v2D*v2D));
}

// This is the exponential factor in the Green's function sum, also
// appearing in the survival prob. and prop. function.
//
// Also here the root is the one refering to [0,a].
Real FirstPassageGreensFunction1DRad::Cn (Real a_n, Real t)
const
{
    const Real D(this->getD());

    return std::exp(-D*a_n*a_n*t);
}

// Calculates the probability of finding the particle inside the domain at
// time t, the survival probability. The domain is from -r to r (r0 is in
// between!!)
Real FirstPassageGreensFunction1DRad::p_survival (Real t) const
{
    const Real D(this->getD());
    const Real v(this->getv());
    const Real vexpo(-v*v*t/4.0/D - v*r0/2.0/D);

    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    if (t == 0.0 || (D == 0.0 && v == 0.0) )
    {
	// if there was no time or no movement the particle was always
	// in the domain
	return 1.0;
    }


    Real an;
    Real sum = 0, term = 0, term_prev = 0;
    int n = 1;

    do
    {
	an = this->a_n(n);
	term_prev = term;
	term = this->Cn(an, t) * this->An(an) * this->Bn(an);
	sum += term;
	n++;
    }
    // TODO: Is 1.0 a good measure for the scale of probability or will this
    // fail at some point?
    while ( fabs(term/sum) > EPSILON*1.0 ||
	fabs(term_prev/sum) > EPSILON*1.0 ||
	n <= MIN_TERMEN);

    return 2.0*exp(vexpo)*sum;
}


// Calculates the probability density of finding the particle at location r at
// time t.
Real FirstPassageGreensFunction1DRad::prob_r (Real r, Real t)
const
{
    const Real a(this->geta());
    const Real D(this->getD());
    const Real v(this->getv());
    const Real h((this->getk()+this->getv()/2.0)/this->getD());
    const Real r0(this->getr0());

    THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    THROW_UNLESS( std::invalid_argument, 0 <= r && r <= a);
    
    const Real vexpo(-v*v*t/D/4.0 + v*(r-r0)/D/2.0);

    // if there was no time or no movement
    if (t == 0 || D == 0)
    {
	// the probability density function is a delta function
	if (r == r0)
	{
	    return INFINITY;
	}
	else
	{
	    return 0.0;
	}
    }

    // if you're looking on the boundary
    if ( fabs (r - a) < EPSILON*a )
    {
	return 0.0;
    }

    Real root_n, root_n_r;
    Real sum = 0, term = 0, prev_term = 0;
    int n=1;

    do
    {
	if ( n >= MAX_TERMEN )
	{
	    std::cerr << "Too many terms needed for GF1DRad::prob_r. N: "
	              << n << std::endl;
	    break;
	}

	root_n = this->a_n(n);
	root_n_r = root_n*r;

	prev_term = term;
	term = Cn(root_n, t) * An(root_n) * (h*sin(root_n_r) + root_n*cos(root_n_r));
	sum += term;

	n++;
    }
    // TODO: PDENS_TYPICAL is now 1e3, is this any good?!
    while (fabs(term/sum) > EPSILON*PDENS_TYPICAL || 
	fabs(prev_term/sum) > EPSILON*PDENS_TYPICAL ||
	n <= MIN_TERMEN );

    return 2.0*exp(vexpo)*sum;
}

// Calculates the probability density of finding the particle at location z at
// timepoint t, given that the particle is still in the domain.
Real
FirstPassageGreensFunction1DRad::calcpcum (Real r, Real t) const
{
    // BEWARE: HERE THERE IS SCALING OF R!
    //const Real r_corr(r/this->l_scale);
    //return prob_r(r_corr, t)/p_survival(t);	// renormalized version, discontinued
    return prob_r(r, t)/p_survival(t);		// p_survival is unscaled!
}

// Calculates the total probability flux leaving the domain at time t
// This is simply the negative of the time derivative of the survival prob.
// at time t [-dS(t')/dt' at t'=t].
Real FirstPassageGreensFunction1DRad::flux_tot (Real t) const
{
    Real an;
    double sum = 0, term = 0, prev_term = 0;
    const Real D(this->getD());
    const Real v(this->getv());
    const Real vexpo(-v*v*t/4.0/D - v*r0/2.0/D);

    int n=1;

    do
    {
	if ( n >= MAX_TERMEN )
	{
	    std::cerr << "Too many terms needed for GF1DRad::flux_tot. N: "
	              << n << std::endl;
	    break;
	}

	an = this->a_n(n);
	prev_term = term;
	term = an * an * Cn(an, t) * this->An(an) * Bn(an);
	n++;
	sum += term;
    }
    while (fabs(term/sum) > EPSILON*PDENS_TYPICAL ||
	fabs(prev_term/sum) > EPSILON*PDENS_TYPICAL ||
	n <= MIN_TERMEN );

    return 2.0*D*exp(vexpo)*sum;
}

// Calculates the probability flux leaving the domain through the radiative
// boundary at time t
Real FirstPassageGreensFunction1DRad::flux_rad (Real t) const
{
    return this->getk()*prob_r(0, t);
}

// Calculates the flux leaving the domain through the radiative boundary as a
// fraction of the total flux. This is the probability that the particle left
// the domain through the radiative boundary instead of the absorbing
// boundary.
Real FirstPassageGreensFunction1DRad::fluxRatioRadTot (Real t) const
{
    return flux_rad(t)/flux_tot(t);
}

// Determine which event has occured, an escape or a reaction. Based on the
// fluxes through the boundaries at the given time. Beware: if t is not a
// first passage time you still get an answer!
EventType
FirstPassageGreensFunction1DRad::drawEventType( Real rnd, Real t )
const
{
    const Real a(this->geta());
    const Real r0(this->getr0());

    THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
    // if t=0 nothing has happened->no event!!
    THROW_UNLESS( std::invalid_argument, t > 0.0 );

    if ( k == 0 || fabs( r0 - a ) < EPSILON*a )
    {
	return IV_ESCAPE;
    }

    const Real fluxratio (this->fluxRatioRadTot(t));

    if (rnd > fluxratio )
    {
	return IV_ESCAPE;
    }
    else
    {
	return IV_REACTION;
    }
}

// This function is needed to cast the math. form of the function
// into the form needed by the GSL root solver.
double FirstPassageGreensFunction1DRad::drawT_f (double t, void *p)
{
    // casts p naar type 'struct drawT_params *'
    struct drawT_params *params = (struct drawT_params *)p;
    Real sum = 0, term = 0, prev_term = 0;
    Real Xn, exponent;
    Real prefactor = params->prefactor;
    int terms = params->terms;

    int n=0;
    do
    {
	if ( n >= terms )
	{
	    std::cerr << "Too many terms needed for GF1DRad::DrawTime. N: "
	              << n << std::endl;
	    break;
	}
	prev_term = term;

	Xn = params->Xn[n];
	exponent = params->exponent[n];
	term = Xn * exp(exponent * t);
	sum += term;
	n++;
    }
    while (fabs(term/sum) > EPSILON*1.0 ||
	fabs(prev_term/sum) > EPSILON*1.0 ||
	n <= MIN_TERMEN );

    // find the intersection with the random number
    return 1.0 - prefactor*sum - params->rnd;
}

// Draws the first passage time from the survival probability,
// using an assistance function drawT_f that casts the math. function
// into the form needed by the GSL root solver.
Real FirstPassageGreensFunction1DRad::drawTime (Real rnd) const
{
    const Real a(this->geta());
    const Real k(this->getk());
    const Real D(this->getD());
    const Real v(this->getv());
    const Real h((this->getk()+this->getv()/2.0)/this->getD());
    const Real r0(this->getr0());

    THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );

    if ( D == 0.0 || a == INFINITY )
    {
	return INFINITY;
    }

    if ( rnd <= EPSILON || a < 0.0 || fabs(r0 - a) < EPSILON*a )
    {
	return 0.0;
    }

    const Real v2D(v/2.0/D);
    const Real exp_av2D(exp(a*v2D));
    // exponent of the prefactor present in case of v<>0; has to be split because it has a t-dep. and t-indep. part
    const Real vexpo_t(-v*v/4.0/D);
    const Real vexpo_pref(-v*r0/2.0/D);

    // the structure to store the numbers to calculate the numbers for 1-S
    struct drawT_params parameters;
    // some temporary variables
    double an = 0;
    double an2, anr0, ana, han;
    double Xn, exponent, prefactor;


    // produce the coefficients and the terms in the exponent and put them
    // in the params structure. This is not very efficient at this point,
    // coefficients should be calculated on demand->TODO
    for (int n=0; n<MAX_TERMEN; n++)
    {
	an = a_n(n+1);	// get the n-th root of tan(root*a)=root/-h (Note: root numbering starts at n=1)
	
	an2 = an * an; 	// an^2
	anr0 = an * r0;	// an * z'
	ana = an * a;	// an * a
	han = h / an;	// h / an
	
	if(v==0)	Xn = (h*sin(anr0) + an*cos(anr0)) / (a*(an2+h*h)+h)
			      * (han + sin(ana) - han*cos(ana)); 
	else		Xn = (h*sin(anr0) + an*cos(anr0)) / (a*(an2+h*h)+h)
			      * (h*k/D - exp_av2D*(an2+h*h)*cos(ana)) / (han * (an2 + v2D*v2D)); 
		  
	exponent = -D*an2 + vexpo_t;

	// store the coefficients in the structure
	parameters.Xn[n] = Xn;
	// also store the values for the exponent
	parameters.exponent[n]=exponent;

    }
    
    // the prefactor of the sum is also different in case of drift<>0 :
    if(v==0)	prefactor = 2.0;
    else	prefactor = 2.0*exp(vexpo_pref);
    parameters.prefactor = prefactor;
    
    // store the random number for the probability
    parameters.rnd = rnd;
    // store the number of terms used
    parameters.terms = MAX_TERMEN;
    parameters.tscale = this->t_scale;

    // Define the function for the rootfinder
    gsl_function F;
    F.function = &FirstPassageGreensFunction1DRad::drawT_f;
    F.params = &parameters;


    // Find a good interval to determine the first passage time in
    // get the distance to absorbing boundary (disregard rad BC)
    const Real dist(a-r0);
    //const Real dist( std::min(r0, a-r0));	// for test purposes
    // construct a guess: MSD = sqrt (2*d*D*t)
    Real t_guess( dist * dist / ( 2.0*D ) );
    // A different guess has to be made in case of nonzero drift to account for the displacement due to it
    // TODO: This does not work properly in this case yet...
    // When drifting towards the closest boundary...
    //if( (r0 >= a/2.0 && v > 0.0) || (r0 <= a/2.0 && v < 0.0) )	t_guess = sqrt(D*D/(v*v*v*v)+dist*dist/(v*v)) - D/(v*v);
    // When drifting away from the closest boundary...
    //if( ( r0 < a/2.0 && v > 0.0) || ( r0 > a/2.0 && v < 0.0) )	t_guess = D/(v*v) - sqrt(D*D/(v*v*v*v)-dist*dist/(v*v));
    
    Real value( GSL_FN_EVAL( &F, t_guess ) );
    Real low( t_guess );
    Real high( t_guess );


    // scale the interval around the guess such that the function straddles
    if( value < 0.0 )
    {
	// if the guess was too low
	do
	{
	    // keep increasing the upper boundary until the
	    // function straddles
	    high *= 10;
	    value = GSL_FN_EVAL( &F, high );

	    if( fabs( high ) >= t_guess * 1e6 )
	    {
		std::cerr << "GF1DRad: Couldn't adjust high. F("
		          << high << ") = " << value << std::endl;
		throw std::exception();
	    }
	}
	while ( value <= 0.0 );
    }
    else
    {
	// if the guess was too high
	// initialize with 2 so the test below survives the first
	// iteration
	Real value_prev( 2 );
	do
	{
	    if( fabs( low ) <= t_guess * 1e-6 ||
	        fabs(value-value_prev) < EPSILON*1.0 )
	    {
		std::cerr << "GF1DRad: Couldn't adjust low. F(" << low << ") = "
		          << value << " t_guess: " << t_guess << " diff: "
		          << (value - value_prev) << " value: " << value
		          << " value_prev: " << value_prev << " rnd: "
		          << rnd << std::endl;
		return low;
	    }
	    value_prev = value;
	    // keep decreasing the lower boundary until the function straddles
	    low *= .1;
	    // get the accompanying value
	    value = GSL_FN_EVAL( &F, low );
	}
	while ( value >= 0.0 );
    }


    // find the intersection on the y-axis between the random number and
    // the function
    // define a new solver type brent
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    // make a new solver instance
    // TODO: incl typecast?
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    const Real t( findRoot( F, solver, low, high, t_scale*EPSILON, EPSILON,
                            "FirstPassageGreensFunction1DRad::drawTime" ) );

    // return the drawn time
    return t;
}

double FirstPassageGreensFunction1DRad::drawR_f (double z, void *p)
{
    // casts p naar type 'struct drawR_params *'
    struct drawR_params *params = (struct drawR_params *)p;
    Real sum = 0, term = 0, prev_term = 0;
    Real an, S_Cn_an, b_an;
    Real v2D = params->H[0];		// = v2D
    Real costerm = params->H[1];	// = k/D
    Real sinterm = params->H[2];	// = h*v2D
    int  terms = params->terms;

    int n = 0;
    do
    {
	if ( n >= terms )
	{
	    std::cerr << "GF1DRad: Too many terms needed for DrawR. N: "
	              << n << std::endl;
	    break;
	}
	prev_term = term;

	S_Cn_an = params->S_Cn_an[n];
	b_an	= params->b_an[n];	// = h/root[n]
	an  = params->an[n];
	term = S_Cn_an * ( costerm - exp(v2D*z)*( costerm*cos(an*z) - (an+sinterm/an)*sin(an*z) ));

	sum += term;
	n++;
    }
    // the function returns a probability (scale is 1, TODO: is this important?)
    while (fabs(term/sum) > EPSILON*1.0 ||
	fabs(prev_term/sum) > EPSILON*1.0 ||
	n <= MIN_TERMEN );

    // Find the intersection with the random number
    return sum - params->rnd;
}

Real
FirstPassageGreensFunction1DRad::drawR (Real rnd, Real t) const
{
    const Real a(this->geta());
    const Real D(this->getD());
    const Real v(this->getv());
    const Real k(this->getk());
    const Real h((this->getk()+this->getv()/2.0)/this->getD());
    const Real r0(this->getr0());

    THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    if (t==0.0 || (D==0.0 && v==0.0) )
    {
	// the trivial case
	//return r0*this->l_scale;	// renormalized version, discontinued
	return r0;
    }
    if ( a<0.0 )
    {
	// if the domain had zero size
	return 0.0;
    }

    // the structure to store the numbers to calculate the numbers for 1-S
    struct drawR_params parameters;
    double an = 0;
    double S_Cn_an;
    double an2, anr0;
    const Real vexpo(-v*v*t/4.0/D - v*r0/2.0/D);	// exponent of the drift-prefactor, same as in survival prob.
    const Real v2D(v/2.0/D);
    const Real S = 2.0*exp(vexpo)/p_survival(t);	// This is a prefactor to every term, so it also contains there
							// exponential drift-prefactor.


    // produce the coefficients and the terms in the exponent and put them
    // in the params structure
    for (int n=0; n<MAX_TERMEN; n++)
    {
	an = a_n(n+1);  // get the n-th root of tan(alfa*a)=alfa/-k
	an2 = an * an; // an^2
	anr0 = an * r0; // an * z'
	S_Cn_an = S * exp(-D*an2*t)
	            * (an*cos(anr0) + h*sin(anr0)) / (a*(an2 + h*h) + h)
		    * an / (an2 + v2D*v2D);

	// store the coefficients in the structure
	parameters.an[n]     = an;
	// also store the values for the exponent
	parameters.S_Cn_an[n]= S_Cn_an;
	parameters.b_an[n]   = h/an;

    }
    // store the random number for the probability
    parameters.rnd = rnd;
    // store the number of terms used
    parameters.terms = MAX_TERMEN;
    
    // also store constant prefactors that appear in the calculation of the
    // r-dependent terms
    parameters.H[0] = v2D;		// appears together with z in one of the prefactors
    parameters.H[1] = k/D;		// further constant terms of the cosine prefactor
    parameters.H[2] = h*v2D;		// further constant terms of the sine prefactor


    // find the intersection on the y-axis between the random number and
    // the function
    gsl_function F;
    F.function = &FirstPassageGreensFunction1DRad::drawR_f;
    F.params = &parameters;

    // define a new solver type brent
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    // make a new solver instance
    // TODO: incl typecast?
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    Real r( findRoot( F, solver, 0.0, a, EPSILON*a, EPSILON,
                            "FirstPassageGreensFunction1DRad::drawR" ) );

    // return the drawn place
    //return r*this->l_scale;	// renormalized version, discontinued
    return r;
}


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
#include "FirstPassageGreensFunction1D.hpp"
#include "Defs.hpp"


// Calculates the probability of finding the particle inside the domain at 
// time t
Real FirstPassageGreensFunction1D::p_survival (Real t) const
{
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    const Real a(this->geta());
    const Real r0(this->getr0());
    const Real D(this->getD());
    const Real v(this->getv());

    if ( a < 0 || fabs (a-r0) < a*EPSILON || r0 < a*EPSILON )
    {
	// The survival probability of a zero domain is zero?
	return 0.0;
    }

    // Set values that are constant in this calculation
    const Real expo(-D*t/(a*a));   // part of the exponent -D n^2 PI^2 t / a^2
    const Real r0_a(r0/a);
    
    // some abbreviations for terms appearing in the sums with drift<>0
    const Real av2D(a*v/2.0/D);
    const Real vexpo(-v*v*t/4.0/D - v*r0/2.0/D);	// exponent of the drift-prefactor
    

    // Initialize summation
    Real sum = 0, term = 0, prev_term = 0;
    Real nPI;


    // Sum
    Real n=1;
    // different calculations depending on whether v=0 or not
    if(v==0.0)	// case without drift (v==0); in this case the summation is simpler, so do the complicated caluclation only if necessary
    {
      do
      {
	  if (n >= MAX_TERMEN )
	  {
	      std::cerr << "Too many terms for p_survival. N: " << n << std::endl;
	      break;
	  }
	  
	  prev_term = term;
	  //nPI = (double)(2*n+1)*M_PI;
	  //term = 2.0 * exp(nPI*nPI*expo) * sin(nPI*r0_a) / nPI;	// TODO: seems not to work in some cases... why?
	  nPI = (double)n*M_PI;
	  term = exp(nPI*nPI*expo) * sin(nPI*r0_a) * (1.0 - cos(nPI)) / nPI;
	  sum += term;
	  n++;
      }
      // Is 1 a good measure or will this fail at some point?
      while (	fabs(term/sum) > EPSILON*1.0 ||
		fabs(prev_term/sum) > EPSILON*1.0 ||
		n < MIN_TERMEN );

      sum = 2.0*sum;	// This is a prefactor of every term, so do only one multiplication here
    }
    else	// case with drift (v<>0)
    {
      do
      {
	  if (n >= MAX_TERMEN )
	  {
	      std::cerr << "Too many terms for p_survival. N: " << n << std::endl;
	      break;
	  }
	  
	  nPI = (double)n*M_PI;
	  prev_term = term;
	  term = exp(nPI*nPI*expo) * (1.0 - cos(nPI)*exp(av2D)) * nPI/(av2D*av2D+nPI*nPI) * sin(nPI*r0_a);
	  sum += term;
	  n++;
      }
      // TODO: Is 1 a good measure or will this fail at some point?
      while (	fabs(term/sum) > EPSILON*1.0 ||
		fabs(prev_term/sum) > EPSILON*1.0 ||
		n < MIN_TERMEN );
      
      sum = 2.0*exp(vexpo) * sum;	// prefactor containing the drift

    }

    return sum;
}

// Calculates the probability density of finding the particle at location r at 
// time t.
Real FirstPassageGreensFunction1D::prob_r (Real r, Real t) const
{
    const Real a(this->geta());
    const Real r0(this->getr0());
    const Real D(this->getD());
    const Real v(this->getv());

    THROW_UNLESS( std::invalid_argument, 0 <= r && r <= a );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

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
    else if ( r < EPSILON || (a-r) < EPSILON || a < 0 )
    {
	return 0.0;
    }

    // Set values that are constant in this calculation
    const Real expo(-D*t/(a*a));
    const Real r_a(r/a);
    const Real r0_a(r0/a);
    const Real vexpo(-v*v*t/4.0/D + v*(r-r0)/2.0/D);	// exponent of the drift-prefactor

    // Initialize summation
    Real nPI;
    Real sum = 0, term = 0, prev_term = 0;

    // Sum
    int n=1;
    do
    {
	if (n >= MAX_TERMEN )
	{
	    std::cerr << "Too many terms for prob_r. N: " << n << std::endl;
	    break;
	}

	prev_term = term;

	nPI = n*M_PI;
	term = exp(nPI*nPI*expo) * sin(nPI*r0_a) * sin(nPI*r_a);
	sum += term;
	n++;
    }
    // TODO: Is 1E3 a good measure for the probability density?!
    while (fabs(term/sum) > EPSILON*PDENS_TYPICAL ||
	fabs(prev_term/sum) > EPSILON*PDENS_TYPICAL ||
	n <= MIN_TERMEN);

    return 2.0/a * exp(vexpo) * sum;
}

// Calculates the probability density of finding the particle at location r at 
// timepoint t, given that the particle is still in the domain.
Real FirstPassageGreensFunction1D::calcpcum (Real r, Real t) const
{
    return prob_r (r, t)/p_survival (t);
}

// Calculates the amount of flux leaving the left boundary at time t
Real FirstPassageGreensFunction1D::leaves(Real t) const
{
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    const Real a(this->geta());
    const Real r0(this->getr0());
    const Real D(this->getD());
    const Real v(this->getv());

    if ( a < 0 || r0 < EPSILON )
    {
	// The flux of a zero domain is zero? Also if the particle 
	// started on the left boundary
	return INFINITY;
    }
    else if ( t < EPSILON*this->t_scale )
    {
	// if t=0.0 the flux must be zero
	return 0.0;
    }


    Real sum = 0, term = 0, prev_term = 0;
    Real nPI;
    const Real D_a_sq(D/(a*a));
    const Real expo(-D_a_sq*t);    // exponent -D n^2 PI^2 t / l^2
    const Real r0_a(r0/a);
    const Real vexpo(-v*v*t/4.0/D - v*r0/2.0/D);
    
    Real n=1;
    do
    {
	if (n >= MAX_TERMEN )
	{
	    std::cerr << "Too many terms for p_survival. N: " << n << std::endl;
	    break;
	}

	nPI = n*M_PI;
	prev_term = term;
	term = nPI * exp(nPI*nPI*expo) * sin(nPI*r0_a);
	sum += term;
	n++;
    }
    // Is PDENS_TYPICAL a good measure or will this fail at some point?
    while (fabs(term/sum) > EPSILON*PDENS_TYPICAL ||
	fabs(prev_term/sum) > EPSILON*PDENS_TYPICAL ||
	n < MIN_TERMEN );

    return 2.0*D_a_sq*exp(vexpo)*sum;
}

// Calculates the amount of flux leaving the right boundary at time t
Real FirstPassageGreensFunction1D::leavea(Real t) const
{
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    const Real a(this->geta());
    const Real r0(this->getr0());
    const Real D(this->getD());
    const Real v(this->getv());

    if ( a < 0 || fabs (a-r0) < EPSILON )
    {
	// The flux of a zero domain is zero? Also if the particle 
	// started on the right boundary.
	return INFINITY;
    }
    else if ( t < EPSILON*this->t_scale )
    {
	// if t=0.0 the flux must be zero
	return 0.0;
    }


    Real sum = 0, term = 0, prev_term = 0;
    Real nPI;
    const Real D_a_sq(D/(a*a));
    const Real expo(-D_a_sq*t);		// exponent -D n^2 PI^2 t / l^2
    const Real r0_a(r0/a);
    const Real vexpo(-v*v*t/4.0/D + v*(a-r0)/2.0/D);
    
    Real n=1;
    do
     {
       if (n >= MAX_TERMEN )
       {
	 std::cerr << "Too many terms for leaves. N: " << n << std::endl;
	 break;
       }
       
       nPI = n*M_PI;
       prev_term = term;
       term = nPI * exp(nPI*nPI*expo) * cos(nPI) * sin(nPI*r0_a);
       sum += term;
       n++;
     }
     // Is PDENS_TYPICAL a good measure or will this fail at some point?
     while (fabs(term/sum) > EPSILON*PDENS_TYPICAL ||
	    fabs(prev_term/sum) > EPSILON*PDENS_TYPICAL ||
	    n < MIN_TERMEN );
     
     return -2.0*D_a_sq*exp(vexpo)*sum;
}

// This draws an eventtype of time t based on the flux through the left (z=0) 
// and right (z=a) boundary. Although not completely accurate, it returns an 
// IV_ESCAPE for an escape through the right boundary and a IV_REACTION for an 
// escape through the left boundary.
EventType
FirstPassageGreensFunction1D::drawEventType( Real rnd, Real t ) const
{
    const Real a(this->geta());
    const Real r0(this->getr0());

    THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
    // if t=0 nothing has happened->no event!!
    THROW_UNLESS( std::invalid_argument, t > 0.0 );

    if ( fabs( r0 - a ) < EPSILON*a )
    {
	// if the particle started on the right boundary
	return IV_ESCAPE;
    }
    else if ( r0 < EPSILON )
    {
	// if the particle started on the left boundary
	return IV_REACTION;
    }

    const Real leaves_s (this->leaves(t));
    const Real leaves_a (this->leavea(t));
    const Real flux_total (leaves_s + leaves_a);
    const Real fluxratio (leaves_s/flux_total);

    if (rnd > fluxratio )
    {
	return IV_ESCAPE;
    }
    else
    {
	return IV_REACTION;
    }
}

// This is a help function that casts the drawT_params parameter structure into
// the right form and calculates the survival probability from it (and returns it).
// The routine drawTime uses this one to sample the next-event time from the
// survival probability using a rootfinder from GSL.
double FirstPassageGreensFunction1D::drawT_f (double t, void *p)
{   
    // casts p to type 'struct drawT_params *'
    struct drawT_params *params = (struct drawT_params *)p;
    Real sum = 0, term = 0, prev_term = 0;
    Real Xn, exponent, prefactor;
    // the maximum number of terms in the params table
    int    terms = params->terms;
    // the timescale used
    Real   tscale = params->tscale;

    int n=0;
    do
    {
	if ( n >= terms )
	{
	    std::cerr << "Too many terms needed for DrawTime. N: "
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
    while (fabs(term/sum) > EPSILON*tscale ||
	fabs(prev_term/sum) > EPSILON*tscale ||
	n <= MIN_TERMEN );

    prefactor = params->prefactor;
    // the intersection with the random number
    return 1.0 - prefactor*sum - params->rnd;
}

// Draws the first passage time from the propensity function.
// Uses the help routine drawT_f and structure drawT_params for some technical
// reasons related to the way to input a function and parameters required by
// the GSL library.
Real FirstPassageGreensFunction1D::drawTime (Real rnd) const
{
    THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );

    const Real a(this->geta());
    const Real r0(this->getr0());
    const Real D(this->getD());
    const Real v(this->getv());
    
    if (D == 0.0 )
    {
	return INFINITY;
    }
    else if ( a < 0.0 || r0 < EPSILON*a || r0 > (1.0 - EPSILON)*a )
    {
	// if the domain had zero size
	return 0.0;
    }

    const Real expo(-D/(a*a));
    const Real r0_a(r0/a);
    // some abbreviations for terms appearing in the sums with drift<>0
    const Real av2D(a*v/2.0/D);
    // exponent of the prefactor present in case of v<>0; has to be split because it has a t-dep. and t-indep. part
    const Real vexpo_t(-v*v/4.0/D);
    const Real vexpo_pref(-v*r0/2.0/D);

    // the structure to store the numbers to calculate the numbers for 1-S
    struct drawT_params parameters;
    Real Xn, exponent, prefactor;
    
    Real nPI;
    
    // Sum
    // construct the coefficients and the terms in the exponent and put them 
    // in the params structure
    int n = 0;
    // a simpler sum has to be computed for the case w/o drift, so distinguish here
    if(v==0)
    {
      do
      {
	  nPI = ((Real)(n+1))*M_PI;	// why n+1 : this loop starts at n=0 (1st index of the arrays), while the sum starts at n=1 !
	  Xn = sin(nPI*r0_a) * (1.0 - cos(nPI)) / nPI; 
	  exponent = nPI*nPI*expo;
	  
	  // store the coefficients in the structure
	  parameters.Xn[n] = Xn;	
	  // also store the values for the exponent
	  parameters.exponent[n]=exponent;
	  n++;
      }
      // modify this later to include a cutoff when changes are small
      while (n<MAX_TERMEN);
    }
    else	// case with drift<>0
    {
     do
      {
	  nPI = ((Real)(n+1))*M_PI;	// why n+1 : this loop starts at n=0 (1st index of the arrays), while the sum starts at n=1 !
	  Xn = (1.0 - cos(nPI)*exp(av2D)) * nPI/(av2D*av2D+nPI*nPI) * sin(nPI*r0_a);
	  exponent = nPI*nPI*expo + vexpo_t;
	  
	  // store the coefficients in the structure
	  parameters.Xn[n] = Xn;	
	  // also store the values for the exponent
	  parameters.exponent[n]=exponent;
	  n++;
      }
      // modify this later to include a cutoff when changes are small
      while (n<MAX_TERMEN);
    }

    // the prefactor of the sum is also different in case of drift<>0 :
    if(v==0)	prefactor = 2.0*exp(vexpo_pref);
    else	prefactor = 2.0;
    parameters.prefactor = prefactor;
    
    // store the random number for the probability
    parameters.rnd = rnd;
    // store the number of terms used
    parameters.terms = MAX_TERMEN;
    parameters.tscale = this->t_scale;

//debugging
/*
for (double t=0; t<0.1; t += 0.0001)
{
    std::cout << t << " " << drawT_f (t, &parameters) << std::endl;
}
*/
    gsl_function F;
    F.function = &drawT_f;
    F.params = &parameters;


    // Find a good interval to determine the first passage time in
    const Real dist( std::min(r0, a-r0));

    // construct a guess: MSD = sqrt (2*d*D*t)
    Real t_guess( dist * dist / ( 2.0 * D ) );
    // A different guess has to be made in case of nonzero drift to account for the displacement due to it
    // When drifting towards the closest boundary...
    if( (r0 >= a/2.0 && v > 0.0) || (r0 <= a/2.0 && v < 0.0) )	t_guess = sqrt(D*D/(v*v*v*v)+dist*dist/(v*v)) - D/(v*v);
    // When drifting away from the closest boundary...
    if( ( r0 < a/2.0 && v > 0.0) || ( r0 > a/2.0 && v < 0.0) )	t_guess = D/(v*v) - sqrt(D*D/(v*v*v*v)-dist*dist/(v*v));
    
    
    Real value( GSL_FN_EVAL( &F, t_guess ) );
    Real low( t_guess );
    Real high( t_guess );

    if( value < 0.0 )
    {
	// scale the interval around the guess such that the function 
	// straddles if the guess was too low
	do
	{
	    // keep increasing the upper boundary until the 
	    // function straddles
	    high *= 10.0;
	    value = GSL_FN_EVAL( &F, high );

	    if( fabs( high ) >= t_guess * 1e6 )
	    {
		std::cerr << "Couldn't adjust high. F(" << high << ") = "
		          << value << std::endl;
		throw std::exception();
	    }
	}
	while ( value <= 0.0 );
    }
    else
    {
	// if the guess was too high initialize with 2 so the test 
	// below survives the first iteration
	Real value_prev( 2.0 );
	do
	{
	    if( fabs( low ) <= t_guess * 1.0e-6 ||
	        fabs(value-value_prev) < EPSILON*this->t_scale )
	    {
		std::cerr << "Couldn't adjust low. F(" << low << ") = "
		          << value << " t_guess: " << t_guess << " diff: "
		          << (value - value_prev) << " value: " << value
		          << " value_prev: " << value_prev << " t_scale: "
		          << this->t_scale << std::endl;
		return low;
	    }

	    value_prev = value;
	    // keep decreasing the lower boundary until the 
	    // function straddles
	    low *= 0.1;
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
    const Real t( findRoot( F, solver, low, high, EPSILON*t_scale, EPSILON,
                            "FirstPassageGreensFunction1D::drawTime" ) );

    // return the drawn time
    return t;
}


// This is a help function that casts the drawR_params parameter structure into
// the right form and calculates the survival probability from it (and returns it).
// The routine drawR uses this function to sample the exit point, making use of the
// GSL root finder to draw the random position.
double FirstPassageGreensFunction1D::drawR_f (double r, void *p)
{   
    struct drawR_params *params = (struct drawR_params *)p;
    double sum = 0, term = 0, prev_term = 0;
    double S_Cn_An, n_l;
    int    terms = params->terms;
    double v2D = params->v2D;	// this is v divided by 2D

    int n=0;
    do
    {
	if (n >= terms )
	{
	    std::cerr << "Too many terms for DrawR. N: " << n << std::endl;
	    break;
	}
	prev_term = term;

	S_Cn_An = params->S_Cn_An[n];
	n_l = params->n_l[n];	// this is n*pi/a
	if(v2D==0.0)	term = S_Cn_An * ( 1.0 - cos(n_l*r) );
	else		term = S_Cn_An * ( 1.0 + exp(v2D*r)*( v2D/n_l*sin(n_l*r) - cos(n_l*r) ));
	  // S_Cn_An contains all expon. prefactors, the 1/S(t) term and all parts
	  // of the terms that do not depend on r.
	  //
	  // In case of zero drift the terms become S_An_Cn * (1 - cos(nPi/a*r))
	  // as it should be.

	sum += term;
	n++;
    }
    // A lengthscale of 1 if implied here; TODO: Is that even true?
    while (fabs(term/sum) > EPSILON ||
	fabs(prev_term/sum) > EPSILON ||
	n <= MIN_TERMEN );

    // find the intersection with the random number
    return sum - params->rnd;
}

// Draws the position of the particle at a given time from p(r,t), assuming 
// that the particle is still in the domain
Real FirstPassageGreensFunction1D::drawR (Real rnd, Real t) const
{
    THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    const Real a(this->geta());
    const Real r0(this->getr0());
    const Real D(this->getD());
    const Real v(this->getv());

    // the trivial case, if there was no movement or the domain was zero
    if ( (D==0.0 && v==0.0) || a<0.0 || t==0.0)
    {
	return r0;
    }
    else
    {
	// if the initial condition is at the boundary, raise an error
	// The particle can only be at the boundary in the ABOVE cases
	THROW_UNLESS( std::invalid_argument,
	              EPSILON <= r0 && r0 <= (a-EPSILON) );
    }
    // else the normal case
    // From here on the problem is well defined


    // structure to store the numbers to calculate numbers for 1-S
    struct drawR_params parameters;
    Real S_Cn_An;
    Real nPI;
    const Real expo (-D*t/(a*a));
    const Real r0_a(r0/a);
    const Real v2D(v/2.0/D);
    const Real vexpo(-v*v*t/4.0/D - v*r0/2.0/D);	// exponent of the drift-prefactor, same as in survival prob.
    const Real S = 2.0*exp(vexpo)/p_survival(t);	// This is a prefactor to every term, so it also contains there
							// exponential drift-prefactor.

    // produce the coefficients and the terms in the exponent and put them 
    // in the params structure
    int n=0;
    do
    {
	nPI = ((Real)(n+1))*M_PI;	    // note: summation starting with n=1, indexing with n=0, therefore we need n+1 here
	
	if(v==0.0)	S_Cn_An = S * exp(nPI*nPI*expo) * sin(nPI*r0_a) / nPI;
	else		S_Cn_An = S * exp(nPI*nPI*expo) * sin(nPI*r0_a) * nPI/(nPI*nPI + a*a*v2D*v2D);
	  // The rest is the z-dependent part, which has to be defined directly in drawR_f(z).
	  // Of course also the summation happens there because the terms now are z-dependent.
	  // The last term originates from the integrated prob. density including drift.
	  //
	  // In case of zero drift this expression becomes: 2.0/p_survival(t) * exp(nPI*nPI*expo) * sin(nPI*r0_a) / nPI
	  
	// also store the values for the exponent, so they don't have to be recalculated in drawR_f
	parameters.S_Cn_An[n]= S_Cn_An;
	parameters.n_l[n]    = nPI/a;
	n++;
    }
    while (n<MAX_TERMEN);

    // store the random number for the probability
    parameters.rnd = rnd ;
    // store the number of terms used
    parameters.terms = MAX_TERMEN;
    
    // store needed constants
    parameters.v2D = v2D;

    // find the intersection on the y-axis between the random number and 
    // the function
    gsl_function F;
    F.function = &drawR_f;
    F.params = &parameters;

    // define a new solver type brent
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    // make a new solver instance
    // TODO: incl typecast?
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    const Real r( findRoot( F, solver, 0.0, a, a*EPSILON, EPSILON,
                            "FirstPassageGreensFunction1D::drawR" ) );

    // return the drawn time
    //return r*(this->l_scale);		// renormalized version, discontinued
    return r;
}


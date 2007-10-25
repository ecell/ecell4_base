
#include <gsl/gsl_math.h>

#include "freeFunctions.hpp"


/**
   Calculates exp( x^2 ) * erfc( x )

   See asymptotic expansion here:
   http://en.wikipedia.org/wiki/Error_function
*/  
const Real expxsq_erfc( const Real x )
{
    Real result;

    const Real xsq( x * x );
    if( x > 26.0 )
    {
	const Real M_1_SQRTPI( M_2_SQRTPI * 0.5 ); 

	const Real x2sq_r( 1.0 / ( 2.0 * xsq ) );  // 2 / (2 x)^2

	/*
	  up to second term in the expansion.
	  abs err ~= 9e-8 at x == 20, 3e-8 at x == 25

	  the third term 
	  - ( 8 / ( x2sq * x2sq * x2sq ) )       
	  and beyond doesn't have a major contribution for large x.
	*/

	result = ( M_1_SQRTPI / x ) * 
	    ( 1.0 - x2sq_r +      // term 1
	      x2sq_r * x2sq_r );  // term 2
    }
    else
    {
	result = exp( xsq ) * erfc( x );
    }

    return result;
}


/**
   W( a, b ) := exp( 2 a b + b^2 ) erfc( a + b )
*/
const Real W( const Real a, const Real b )
{
    // exp( 2 a b + b^2 ) erfc( a + b ) == 
    //               exp( - a^2 ) exp( ( a + b )^2 ) erfc( a + b )
    return exp( - a * a ) * expxsq_erfc( a + b );
}

const Real 
__p_irr( const Real r,
         const Real t,
         const Real r0,
         const Real kf,
         const Real D,
         const Real sigma,
         const Real alpha )
{
    //  printf("irrp %g %g %g\n",r,r0,t);
    const Real sqrtD( sqrt( D ) );

    const Real Dt4( 4.0 * D * t );
    const Real r_plus_r0_minus_2sigma( r + r0 - 2.0 * sigma );

    const Real num1( exp( - gsl_pow_2( r - r0 ) / Dt4 ) );
    const Real num2( exp( - gsl_pow_2( r_plus_r0_minus_2sigma ) / Dt4 ) );
    const Real num3( W( r_plus_r0_minus_2sigma / sqrt( Dt4 ), 
			alpha * sqrt( t ) ) );

    const Real num( ( num1 + num2 ) / sqrt( 4.0 * M_PI * t ) -  alpha * num3 );

    const Real den( 4.0 * M_PI * r * r0 * sqrtD );

    const Real result( num / den );

    const Real jacobian( 4.0 * M_PI * r * r );

    return result * jacobian;
}

const Real 
p_irr( const Real r,
        const Real t,
        const Real r0,
        const Real kf,
        const Real D,
        const Real sigma )
{
    const Real kD( 4.0 * M_PI * sigma * D );
    const Real alpha( ( 1.0 + ( kf / kD ) ) * ( sqrt( D ) / sigma ) );

    const Real p( __p_irr( r, t, r0, kf, D, sigma, alpha ) );

    return p;
}


const Real S_irr( const Real t, const Real r0,
                  const Real kf, const Real D, const Real sigma )
{
    const Real kD( 4.0 * M_PI * sigma * D );
    const Real alpha( ( 1.0 + ( kf / kD ) ) * ( sqrt( D ) / sigma ) );

    const Real sqrtt( sqrt( t ) );

    const Real p( __p_reaction_irr( sqrtt, r0, kf, D, sigma, alpha, kD ) );

    return 1.0 - p;
}

const Real 
__p_reaction_irr( const Real sqrtt, const Real r0,
                  const Real kf, const Real D, const Real sigma,
                  const Real alpha, const Real kD )
{
    const Real sqrtD( sqrt( D ) );

    const Real r0_m_sigma_over_sqrt4D_t( ( r0 - sigma ) 
					 / ( ( sqrtD + sqrtD ) * sqrtt ) );

    const Real Wf( W( r0_m_sigma_over_sqrt4D_t, alpha * sqrtt ) );
    const Real factor( sigma * kf / ( r0 * ( kf + kD ) ) );

    return factor * ( erfc( r0_m_sigma_over_sqrt4D_t ) - Wf );
}

/*
const Real S_irr_deriv( const Real tsqrt, 
                        const Real r0 ) const
{
    const Real Sigma( this->getSigma() );
    const Real D( this->getD() );
    const Real alpha( this->getalpha() );
    const Real kD( this->getkD() );
    const Real kf( this->getkf() );

    const Real sqrtD( sqrt( D ) );
    const Real sqrtPI( sqrt( M_PI ) );

    const Real r0_m_Sigma_t_over_sqrt4D( ( r0 - Sigma ) * tsqrt / 
					 ( sqrtD + sqrtD ) );
    const Real Wf( W( r0_m_Sigma_t_over_sqrt4D, alpha * tsqrt ) );

    const Real num1( sqrtD * exp( - gsl_pow_2( r0_m_Sigma_t_over_sqrt4D ) ) );
    const Real num2( ( sqrtPI * tsqrt * ( alpha * sqrtD + r0 - Sigma ) ) * Wf );

    const Real factor( ( alpha + alpha ) * kf * Sigma /
		       ( sqrtPI * sqrtD * r0 * ( kf + kD ) ) );
  
    return ( num1 - num2 ) * factor;
}

void
S_irr_fdf( const Real tsqrt, 
					 const Real r0,
					 Real* const f, Real* const df ) const
{
    const Real kD( this->getkD() );
    const Real kf( this->getkf() );
    const Real Sigma( this->getSigma() );
    const Real D( this->getD() );
    const Real alpha( this->getalpha() );

    const Real sqrtD( sqrt ( D ) );

    const Real r0_m_Sigma_over_sqrt4D( ( r0 - Sigma ) / ( sqrtD + sqrtD ) );
    const Real factor( Sigma * kf / ( r0 * ( kf + kD ) ) );

    {
	const Real r0_m_Sigma_over_sqrt4D_t( r0_m_Sigma_over_sqrt4D / tsqrt );
	const Real Wf( W( r0_m_Sigma_over_sqrt4D_t, alpha * tsqrt ) );

	*f = factor * ( erfc( r0_m_Sigma_over_sqrt4D_t ) - Wf );
    }

    {
	const Real r0_m_Sigma_t_over_sqrt4D( r0_m_Sigma_over_sqrt4D * tsqrt );
	const Real Wdf( W( r0_m_Sigma_t_over_sqrt4D, alpha * tsqrt ) );
	const Real sqrtPI( sqrt( M_PI ) );

	const Real dfnum1( sqrtD * 
			   exp( - gsl_pow_2( r0_m_Sigma_t_over_sqrt4D ) ) );
	const Real dfnum2( ( sqrtPI * tsqrt * ( alpha * sqrtD + r0 - Sigma ) ) 
			   * Wdf );
    
	const Real dffactor( ( alpha * M_2_SQRTPI / sqrtD ) * factor );
    
	*df = ( dfnum1 - dfnum2 ) * dffactor;
    }
}
*/

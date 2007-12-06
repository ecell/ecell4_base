
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


const Real p_theta_free( const Real theta, const Real r, const Real r0, 
                         const Real t, const Real D )
{
    Real sin_theta;
    Real cos_theta;
    sincos( theta, &sin_theta, &cos_theta );

    const Real Dt4( 4.0 * D * t );
    const Real Dt4Pi( Dt4 * M_PI );

    const Real term1( exp( - ( r * r - 2.0 * cos_theta * r * r0 + r0 * r0 ) / 
                           Dt4 ) );
    const Real term2( 1.0 / sqrt( Dt4Pi * Dt4Pi * Dt4Pi ) );

    return term1 * term2 * sin_theta; // jacobian
}


const Real g_bd( const Real r0, const Real sigma, const Real t, const Real D )
{
    const Real Dt4( 4.0 * D * t );
    const Real mDt4_r( - 1.0 / Dt4 );
    const Real sqrtDt4( sqrt( Dt4 ) );
    const Real sqrtDt4_r( 1.0 / sqrtDt4 );
    const Real sqrtPi( sqrt( M_PI ) );

    const Real r0ps( r0 + sigma );
    const Real r0ms( r0 - sigma );

    const Real term1( ( exp( r0ps * r0ps * mDt4_r ) - 
                        exp( r0ms * r0ms * mDt4_r ) ) * sqrtDt4 / 
                      ( sqrtPi * r0 ) );
    const Real term2( erf( r0ps * sqrtDt4_r ) - erf( r0ms * sqrtDt4_r ) );

    return 0.5 * ( term1 + term2 ) * 4.0 * M_PI * r0 * r0;
}
    
const Real I_bd( const Real sigma, const Real t, const Real D )
{
    const Real sqrtPi( sqrt( M_PI ) );

    const Real Dt( D * t );
    const Real Dt2( 2.0 * Dt );
    const Real sqrtDt( sqrt( Dt ) );
    const Real sigmasq( sigma * sigma );
    const Real sigmasq_over_Dt( sigmasq / Dt );

    const Real exp_sigmasq_over_Dt( exp( - sigmasq_over_Dt ) );

    const Real term1( 4.0 / 3.0 * sqrtPi );
    const Real term2( sigmasq - Dt2 );
    const Real term3( - 3.0 * sigmasq + Dt2 );
    const Real term4( sqrtPi * sigmasq * sigma * erfc( sigma / sqrtDt ) );

    const Real result( term1 * ( ( - sqrtDt *
                                   ( term2 * exp_sigmasq_over_Dt + term3 ) )
                                 + term4 ) );
    
    return result;
}

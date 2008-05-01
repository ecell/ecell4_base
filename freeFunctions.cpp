#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

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

    const Real p( __p_reaction_irr( t, r0, kf, D, sigma, alpha, kD ) );

    return 1.0 - p;
}

const Real 
__p_reaction_irr( const Real t, const Real r0,
                  const Real kf, const Real D, const Real sigma,
                  const Real alpha, const Real kD )
{
    const Real sqrtt( sqrt( t ) );
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

const Real ip_theta_free( const Real theta, const Real r, const Real r0,
                          const Real t, const Real D )
{
    const Real Dt( D * t );
    const Real Dt2( Dt + Dt );
    const Real rr0( r * r0 );

    const Real rr0_over_2Dt( rr0 / Dt2 );

    const Real rsqr0sq_over_4Dt( ( r * r + r0 * r0 ) / ( Dt2 + Dt2 ) );

    const Real term1( expm1( rr0_over_2Dt 
                             - rsqr0sq_over_4Dt ) );
    const Real term2( expm1( rr0_over_2Dt * cos( theta ) 
                             - rsqr0sq_over_4Dt ) );

    const Real den( 4.0 * sqrt( M_PI * M_PI * M_PI * Dt ) * rr0 );

    return ( term1 - term2 ) / den;
}



const Real g_bd( const Real r, const Real sigma, const Real t, const Real D )
{
    const Real Dt4( 4.0 * D * t );
    const Real mDt4_r( - 1.0 / Dt4 );
    const Real sqrtDt4( sqrt( Dt4 ) );
    const Real sqrtDt4_r( 1.0 / sqrtDt4 );
    const Real sqrtPi( sqrt( M_PI ) );

    const Real rps( r + sigma );
    const Real rms( r - sigma );

    const Real term1( ( exp( rps * rps * mDt4_r ) - 
                        exp( rms * rms * mDt4_r ) ) * sqrtDt4 / 
                      ( sqrtPi * r ) );
    const Real term2( erf( rps * sqrtDt4_r ) - erf( rms * sqrtDt4_r ) );

    return 0.5 * ( term1 + term2 ) * r * r;
}
    
const Real I_bd( const Real sigma, const Real t, const Real D )
{
    const Real sqrtPi( sqrt( M_PI ) );

    const Real Dt( D * t );
    const Real Dt2( Dt + Dt );
    const Real sqrtDt( sqrt( Dt ) );
    const Real sigmasq( sigma * sigma );

    const Real term1( 1.0 / ( 3.0 * sqrtPi ) );
    const Real term2( sigmasq - Dt2 );
    const Real term3( Dt2 - 3.0 * sigmasq );
    const Real term4( sqrtPi * sigmasq * sigma * erfc( sigma / sqrtDt ) );

    const Real result( term1 * ( - sqrtDt *
                                 ( term2 * exp( - sigmasq / Dt ) + term3 )
                                 + term4 ) );
    
    return result;
}


const Real I_bd_r( const Real r, const Real sigma, const Real t, const Real D )
{
    const Real sqrtPi( sqrt( M_PI ) );

    const Real Dt( D * t );
    const Real Dt2( Dt + Dt );
    const Real Dt4( Dt2 + Dt2 );
    const Real sqrtDt( sqrt( Dt ) );
    const Real sqrtDt4( sqrt( Dt4 ) );
    const Real sigmasq( sigma * sigma );

    const Real sigmacb( sigmasq * sigma );
    const Real rcb( gsl_pow_3( r ) );

    const Real rsigma( r * sigma );

    const Real rps_sq( gsl_pow_2( r + sigma ) );
    const Real rms_sq( gsl_pow_2( r - sigma ) );

    const Real term1( - 2.0 * sqrtDt / sqrtPi );
    const Real term2( exp( - sigmasq / Dt ) * ( sigmasq - Dt2 ) );
    const Real term3( - exp( - rps_sq / Dt4 ) * ( rms_sq + rsigma - Dt2 ) );
    const Real term4( exp( - rms_sq / Dt4 ) * ( rps_sq - rsigma - Dt2 ) );
    const Real term5( - sigmasq * 3.0 + Dt2 );

    const Real term6( ( sigmacb - rcb ) * erf( ( r - sigma ) / sqrtDt4 ) );
    const Real term7( - ( sigmacb + sigmacb ) * erf( sigma / sqrtDt ) );
    const Real term8( ( sigmacb + rcb ) * erf( ( r + sigma ) / sqrtDt4 ) );

    /* FIXME: expm1, erfc?
    printf("%g %g %g %g %g %g\n", 
           expm1( - sigmasq / Dt ),
           expm1( - rps_sq / Dt4 ),
           expm1( - rms_sq / Dt4 ),
           erf( ( r - sigma ) / sqrtDt4 ),
           erf( sigma / sqrtDt ),
           erf( ( r + sigma ) / sqrtDt4 ) );*/


    const Real result( ( term1 * ( term2 + term3 + term4 + term5 )
                         // + sigmasq + rsigma + rsigma - Dt2 )//expm1
                         + term6 + term7 + term8 ) / 6.0 );
    
    return result;
}


struct g_bd_params
{ 
    const Real sigma;
    const Real t;
    const Real D;
    const Real target;
};


static const Real I_gbd_r_F( const Real r,
                             const g_bd_params* const params )
{
    const Real sigma( params->sigma );
    const Real t( params->t );
    const Real D( params->sigma );
    const Real target( params->target );

    //printf("I %g\n",I_bd_r( r, sigma, t, D ) - target);
    return I_bd_r( r, sigma, t, D ) - target;
}

const Real drawR_gbd( const Real rnd, const Real sigma, 
                      const Real t, const Real D )
{
    const Real I( I_bd( sigma, t, D ) );
//    printf("II %g\n", I );
    g_bd_params params = { sigma, t, D, rnd * I };

    gsl_function F =
    {
        reinterpret_cast<typeof(F.function)>( &I_gbd_r_F ),
        &params
    };

    Real low( sigma );
    Real high( sigma + 10.0 * sqrt ( 6.0 * D * t ) );

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, low, high );

    const unsigned int maxIter( 100 );

    unsigned int i( 0 );
    while( true )
    {
	gsl_root_fsolver_iterate( solver );

	low = gsl_root_fsolver_x_lower( solver );
	high = gsl_root_fsolver_x_upper( solver );
	int status( gsl_root_test_interval( low, high, 1e-18, 1e-12 ) );

	if( status == GSL_CONTINUE )
	{
	    if( i >= maxIter )
	    {
		gsl_root_fsolver_free( solver );
		std::cerr << "drawR_gbd: failed to converge." << std::endl;
		throw std::exception();
	    }
	}
	else
	{
	    break;
	}

	++i;
    }
  
    gsl_root_fsolver_free( solver );

    return low;
}


#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>


#include "FreePairGreensFunction.hpp"




const Real 
FreePairGreensFunction::p_r( const Real r, 
                             const Real r0, 
                             const Real t ) const
{
    const Real D( getD() );
    const Real Dt( D * t );
    const Real Dt4( 4.0 * Dt );
    const Real rr04( 4.0 * r * r0 );

    const Real mrr0sq_over_4Dt( - gsl_pow_2( r + r0 ) / Dt4 );

    const Real num1( expm1( mrr0sq_over_4Dt ) );
    const Real num2( expm1( mrr0sq_over_4Dt + rr04 / Dt4 ) );

    const Real den( rr04 * sqrt( M_PI * M_PI * M_PI * Dt ) );

    const Real jacobian( 2.0 * r * r * M_PI );

    return jacobian * ( - num1 + num2 ) / den;
}

const Real 
FreePairGreensFunction::ip_r( const Real r, 
                              const Real r0, 
                              const Real t ) const
{
    const Real D( getD() );
    const Real Dt4( 4.0 * D * t );
    const Real Dt4r( 1.0 / Dt4 );
    const Real Dt4sqrt( sqrt( Dt4 ) );
    const Real Dt4sqrtr( 1.0 / Dt4sqrt );

    const Real num1a( exp( - gsl_pow_2( r - r0 ) * Dt4r ) );
    const Real num1b( exp( - gsl_pow_2( r + r0 ) * Dt4r ) );
    const Real den1( r0 * sqrt( M_PI ) );

    const Real term1( Dt4sqrt * ( - num1a + num1b ) / den1 );

    const Real term2( erf( ( r - r0 ) * Dt4sqrtr ) );
    const Real term3( erf( ( r + r0 ) * Dt4sqrtr ) );

    return ( term1 + term2 + term3 ) * .5;
}
    
const Real 
FreePairGreensFunction::p_theta( const Real theta,
                                 const Real r,
                                 const Real r0,
                                 const Real t ) const
{
    const Real D( getD() );
    const Real Dt4( 4.0 * D * t );
    const Real Dt4PI( Dt4 * M_PI );

    Real sin_theta;
    Real cos_theta;
    sincos( theta, &sin_theta, &cos_theta );
    const Real p( ( 1.0 / std::sqrt( gsl_pow_3( Dt4PI ) ) ) *
                  exp( - ( r * r + r0 * r0 - 2.0 * r * r0 * cos_theta ) 
                       / Dt4 ) );

    return p * sin_theta * 0.5;
}

const Real 
FreePairGreensFunction::ip_theta( const Real theta,
                                  const Real r,
                                  const Real r0,
                                  const Real t ) const
{
    const Real Dt( getD() * t );
    const Real Dt2( Dt + Dt );
    const Real rr0( r * r0 );

    const Real rr0_over_2Dt( rr0 / ( Dt2 ) );

    const Real rsqr0sq_over_4Dt( ( r * r + r0 * r0 ) / ( Dt2 + Dt2 ) );

    const Real term1( expm1( rr0_over_2Dt 
                             - rsqr0sq_over_4Dt ) );
    const Real term2( expm1( rr0_over_2Dt * cos( theta ) 
                             - rsqr0sq_over_4Dt ) );

    const Real den( 2.0 * 4.0 * sqrt( M_PI * M_PI * M_PI * Dt ) * rr0 );

    return ( term1 - term2 ) / den;
}

const Real
FreePairGreensFunction::ip_r_F( const Real r,
                                const ip_r_params* params )
{
    const FreePairGreensFunction* const gf( params->gf ); 
    const Real r0( params->r0 );
    const Real t( params->t );
    const Real value( params->value );

    return gf->ip_r( r, r0, t ) - value;
}



const Real 
FreePairGreensFunction::drawR( const Real rnd, 
                               const Real r0, 
                               const Real t ) const
{
    // input parameter range checks.
    THROW_UNLESS( std::invalid_argument, rnd <= 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r0 >= 0.0 );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    // t == 0 means no move.
    if( t == 0.0 )
    {
	return r0;
    }

    ip_r_params params = { this, r0, t, rnd };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &ip_r_F ),
	    &params 
	};

    const Real max_r( H * sqrt( 6.0 * getD() * t ) + r0 );

    if( GSL_FN_EVAL( &F, max_r ) < 0.0 )
    {
        return max_r;
    }

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, 0.0, max_r );

    const unsigned int maxIter( 100 );

    unsigned int i( 0 );
    while( true )
    {
	gsl_root_fsolver_iterate( solver );
	const Real low( gsl_root_fsolver_x_lower( solver ) );
	const Real high( gsl_root_fsolver_x_upper( solver ) );
	const int status( gsl_root_test_interval( low, high, 1e-15, 
						  this->TOLERANCE ) );

	if( status == GSL_CONTINUE )
	{
	    if( i >= maxIter )
	    {
		gsl_root_fsolver_free( solver );
		std::cerr << "drawR: failed to converge." << std::endl;
		throw std::exception();
	    }
	}
	else
	{
	    break;
	}

	++i;
    }
  
    //printf("%d\n", i );

    const Real r( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );
    
    return r;
}

const Real
FreePairGreensFunction::ip_theta_F( const Real theta,
                                    const ip_theta_params* params )
{
    const FreePairGreensFunction* const gf( params->gf ); 
    const Real r( params->r );
    const Real r0( params->r0 );
    const Real t( params->t );
    const Real value( params->value );

    return gf->ip_theta( theta, r, r0, t ) - value;
}



const Real 
FreePairGreensFunction::drawTheta( const Real rnd,
                                   const Real r, 
                                   const Real r0, 
                                   const Real t ) const
{
    // input parameter range checks.
    THROW_UNLESS( std::invalid_argument, rnd <= 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r0 >= 0.0 );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    // t == 0 means no move.
    if( t == 0.0 )
    {
	return 0.0;
    }

    const Real ip_theta_pi( ip_theta( M_PI, r, r0, t ) );

    ip_theta_params params = { this, r, r0, t, rnd * ip_theta_pi };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &ip_theta_F ),
	    &params 
	};

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, 0.0, M_PI );

    const unsigned int maxIter( 100 );

    unsigned int i( 0 );
    while( true )
    {
	gsl_root_fsolver_iterate( solver );
	const Real low( gsl_root_fsolver_x_lower( solver ) );
	const Real high( gsl_root_fsolver_x_upper( solver ) );
	const int status( gsl_root_test_interval( low, high, 1e-15, 
						  this->TOLERANCE ) );

	if( status == GSL_CONTINUE )
	{
	    if( i >= maxIter )
	    {
		gsl_root_fsolver_free( solver );
		std::cerr << "drawTheta: failed to converge." << std::endl;
		throw std::exception();
	    }
	}
	else
	{
	    break;
	}

	++i;
    }
  
    //printf("%d\n", i );

    const Real theta( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );
    

    return theta;
}

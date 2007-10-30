#include <iostream>
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


#include "FirstPassageGreensFunction.hpp"




/**
  EllipticTheta[4,0,q]

  Efficiently calculate EllipticTheta[4,0,q] for q < 1.0.
*/

static const Real ellipticTheta4Zero( const Real q )
{
    THROW_UNLESS( std::invalid_argument, fabs( q ) <= 1.0 );
    
    // et4z( 1 - 1e4 ) ~= 7.2e-23
    // et4z( 1e-15 ) ~= 1 - 2e-15
    // et4z( 1e-16 ) ~= 1 - 2.2e-16
    // et4z( 1e-17 ) ~= 1 - (zero)

    const Integer N( 100 );
    Real value( 1.0 );

    Real q_n( q );
    Real q_2n( 1.0 );

    for( Integer n( 1 ); n <= N; ++n )
    {
	const Real term2( 1.0 - q_2n * q );  // q^(2n-1) = (q^(n-1))^2 * q

	q_2n = q_n * q_n;

	const Real term1( 1.0 - q_2n ); // q^2n

	const Real term( term1 * term2 * term2 );
	const Real value_prev( value );
	value *= term;
      
	// here only absolute error is checked because it is good enough
	// for our use in p_survival().
	if( fabs( value - value_prev ) < 1e-8 ) 
	{
	    // normal exit.
	    return value;
	}

	q_n *= q;  // q_(++n)
    }

    std::cerr << "WARNING: ellipticTheta4Zero: didn't converge." << std::endl;
    return value;
}


const Real 
FirstPassageGreensFunction::p_survival( const Real t ) const
{
    const Real D( getD() );
    const Real a( geta() );
    const Real asq( a * a );
    const Real PIsq( M_PI * M_PI );

    const Real q( - D * PIsq * t / asq );

    return 1.0 - ellipticTheta4Zero( exp( q ) );
} 




const Real 
FirstPassageGreensFunction::p_free_int( const Real r, const Real t ) const
{
    const Real D( getD() );
    const Real Dt( D * t );
    const Real sqrtD( sqrt( D ) );
    const Real sqrtt( sqrt( t ) );
    const Real sqrtPI( sqrt( M_PI ) );
    const Real sqrtPI2Dt15( 2.0 * sqrtPI * sqrt( Dt * Dt * Dt ) );

    return ( sqrtPI2Dt15 * erf( r / ( 2.0 * sqrtD * sqrtt ) )
	     - 2.0 * Dt * r * exp( - r * r / ( 4.0 * Dt ) ) ) / sqrtPI2Dt15;
}

const Real 
FirstPassageGreensFunction::p_r_int( const Real r, 
                                     const Real t ) const
{
    Real value( 0.0 );

    const Real a( geta() );
    const Real p_free( this->p_free_int( r, t ) );

    // p_r_int is always smaller than p_free.
    if( fabs( p_free ) < CUTOFF )
    {
	return 0.0;
    }

    const Real D( getD() );
    const Real asq( a * a );
    const Real PIsq( M_PI * M_PI );

    const Real PIr( M_PI * r );
    const Real PIr_a( PIr / a );
    const Real DtPIsq_asq( D * t * PIsq / asq );

    const Real factor( 2.0 / ( a * M_PI ) );

    const Integer N( 10000 );
    long int n( 1 );
    while( true )
    {
	const Real term1( exp( - n * n * DtPIsq_asq ) );
      
	const Real angle_n( n * PIr_a );
	Real sin_n;
	Real cos_n;
	sincos( angle_n, &sin_n, &cos_n );
	const Real term2( a * sin_n );
	const Real term3( n * PIr * cos_n );

	const Real term( term1 * ( term2 - term3 ) / n );
	value += term;

	if( fabs( value ) * CUTOFF >= fabs( term ) )
	{
	    break;
	}

	if( n > N )
	{
	    std::cerr << "p_r_int: didn't converge; " << n << " " << value 
		      << std::endl;
	    break;
	}

	++n;
    }

    //  printf( "value: %g, Dt/a^2 %g\tfree %g\n", value*factor, D*t/(a*a),p_free_int( r, t ) );
    return value * factor;
} 





const Real 
FirstPassageGreensFunction::p_r_fourier( const Real r, const Real t ) const 
{
    Real value( 0.0 );

    const Real D( getD() );
    const Real a( geta() );
    const Real asq( a * a );
    const Real PIsq( M_PI * M_PI );

    const Integer N( 100 );

    long int n( 1 );
    while( true )
    {
	const Real term1( exp( - ( PIsq * r * r + asq * n*n ) / 
			       ( 4.0 * D * PIsq * t ) ) );

	const Real term2( M_PI * r * 
			  exp( gsl_sf_lncosh( a * r * n / 
					      ( 2.0 * D * M_PI * t ) ) ) );

	const Real term3( a * n *
			  exp( gsl_sf_lnsinh( a * r * n / 
					      ( 2.0 * D * M_PI * t ) ) ) );


	const Real term( term1 * r * ( term2 - term3 ) );
	value += term;

	//      printf("%d %g %g %g %g\n", n, value, term, term2, term3 );

	if( fabs( value ) * 1e-8 > fabs( term ) )
	{
	    break;
	}

	if( n > N )
	{
	    std::cerr << "p_r_fourier: didn't converge; " << n << " " << value 
		      << std::endl;
	    break;
	}

	++n;
    }

    const Real factor( 1.0 / ( sqrt( 2 ) * PIsq * pow(D * t, 1.5) ) );

    return value * factor;
} 


const Real
FirstPassageGreensFunction::p_survival_F( const Real t,
					  const p_survival_params* params )
{
    const FirstPassageGreensFunction* const gf( params->gf ); 
    const Real rnd( params->rnd );

    return rnd - gf->p_survival( t );
}



const Real 
FirstPassageGreensFunction::drawTime( const Real rnd ) const
{
    THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );

    const Real a( geta() );

    if( a == 0.0 )
    {
	return 0.0;
    }

    if( getD() == 0.0 )
    {
        return INFINITY;
    }


    p_survival_params params = { this, rnd };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &p_survival_F ),
	    &params 
	};

    Real low( 1e-6 );
    Real high( 1.0 );

    // adjust low to make sure that f( low ) and f( high ) straddle.
    while( GSL_FN_EVAL( &F, high ) <= 0.0 )
    {
	//printf("drawTime: adjusting high: %g\n",high);
	high *= 10;
	if( fabs( high ) >= INFINITY )
	{
            return INFINITY;
	}
    }
    while( GSL_FN_EVAL( &F, low ) >= 0.0 )
    {
	//printf("drawTime: adjusting low: %g\n",low);
	low *= .1;
	if( fabs( low ) <= 1e-50 )
	{
	    std::cerr << "Couldn't adjust low. (" << low <<
		      ")" << std::endl;
	    throw std::exception();
	    
	}
    }


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
		std::cerr << "drawTime: failed to converge." << std::endl;
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


const Real
FirstPassageGreensFunction::p_r_F( const Real r,
				   const p_r_params* params )
{
    const FirstPassageGreensFunction* const gf( params->gf ); 
    const Real t( params->t );
    const Real St( params->St );
    const Real rnd( params->rnd );

    return gf->p_r_int( r, t ) - rnd * St;
}


const Real 
FirstPassageGreensFunction::drawR( const Real rnd, const Real t ) const 
{
    THROW_UNLESS( std::invalid_argument, rnd <= 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    const Real a( geta() );

    if( a == 0.0 || t == 0.0 || getD() == 0.0 )
    {
        return 0.0;
    }

    const Real psurv( p_survival( t ) ); 

    if( psurv == 0.0 )
    {
	printf("p_survival = 0.0\n");
    }

    p_r_params params = { this, t, psurv, rnd };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &p_r_F ),
	    &params 
	};

    Real low( 0.0 );
    Real high( a );

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
  
    gsl_root_fsolver_free( solver );

    return low;
}




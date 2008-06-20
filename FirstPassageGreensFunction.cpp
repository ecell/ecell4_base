#include <sstream>
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

#include "findRoot.hpp"

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

    const Integer N( 1000 );
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
	// for our use.  (it's compared with 1 in p_survival).
	if( fabs( value - value_prev ) < 1e-23 ) 
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
FirstPassageGreensFunction::p_int_r_free( const Real r, const Real t ) const
{
    const Real D( getD() );
    const Real Dt( D * t );
    const Real sqrtDt( sqrt( Dt ) );
    const Real sqrtPI( sqrt( M_PI ) );

    return erf( r / ( sqrtDt + sqrtDt ) )
        - r * exp( - r * r / ( 4.0 * Dt ) ) / ( sqrtPI * sqrtDt );
}

const Real 
FirstPassageGreensFunction::p_int_r( const Real r, 
                                     const Real t ) const
{
    Real value( 0.0 );

    const Real a( geta() );
    const Real p_free( this->p_int_r_free( r, t ) );

    // p_int_r is always smaller than p_free.
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

    const Real maxn( ( a / M_PI ) * sqrt( log( exp( DtPIsq_asq ) / CUTOFF ) / 
                                          ( D * t ) ) );

    const Integer N( std::min( static_cast<Integer>( ceil( maxn ) + 1 ),
                               10000 ) );

    for( int n( 1 ); n <= N; ++n )
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

    if( getD() == 0.0 || geta() == INFINITY )
    {
        return INFINITY;
    }

    if( a == 0.0 )
    {
	return 0.0;
    }

    p_survival_params params = { this, rnd };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &p_survival_F ),
	    &params 
	};

    const Real t_guess( a * a / ( 6. * D ) );

    Real low( t_guess );
    Real high( t_guess );

    const Real value( GSL_FN_EVAL( &F, t_guess ) );

    if( value < 0.0 )
    {
        high *= 10;

        while( 1 )
        {
            const Real high_value( GSL_FN_EVAL( &F, high ) );
            
            if( high_value >= 0.0 )
            {
                break;
            }

            if( fabs( high ) >= t_guess * 1e6 )
            {
                std::cerr << "Couldn't adjust high. F(" << high <<
                    ") = " << GSL_FN_EVAL( &F, high ) << "; " <<
                    ", " << dump() << std::endl;
                throw std::exception();
            }
            high *= 10;
        }
    }
    else
    {
        Real low_value_prev( value );
        low *= .1;

        while( 1 )
        {
            const Real low_value( GSL_FN_EVAL( &F, low ) );
            
            if( low_value <= 0.0 )
            {
                break;
            }
            
            if( fabs( low ) <= t_guess * 1e-6 ||
                fabs( low_value - low_value_prev ) < CUTOFF )
            {
                std::cerr << "Couldn't adjust low.  Returning low (= "
                          << low << "); F(" << low <<
                    ") = " << GSL_FN_EVAL( &F, low )
                          << dump() << std::endl;
                return low;
            }
            low_value_prev = low_value;
            low *= .1;
        }
    }


    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    const Real t( findRoot( F, solver, low, high, 1e-18, 1e-12,
                            "FirstPassageGreensFunction::drawTime" ) );

    gsl_root_fsolver_free( solver );

    return t;
}

const Real
FirstPassageGreensFunction::p_r_free_F( const Real r,
                                        const p_r_params* params )
{
    const FirstPassageGreensFunction* const gf( params->gf ); 
    const Real t( params->t );
    const Real target( params->target );

    return gf->p_int_r_free( r, t ) - target;
}


const Real
FirstPassageGreensFunction::p_r_F( const Real r,
				   const p_r_params* params )
{
    const FirstPassageGreensFunction* const gf( params->gf ); 
    const Real t( params->t );
    const Real target( params->target );

    return gf->p_int_r( r, t ) - target;
}



const Real 
FirstPassageGreensFunction::drawR( const Real rnd, const Real t ) const 
{
    THROW_UNLESS( std::invalid_argument, rnd <= 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    const Real a( geta() );
    const Real D( getD() );

    if( a == 0.0 || t == 0.0 || D == 0.0 )
    {
        return 0.0;
    }

    const Real psurv( p_survival( t ) ); 
    //const Real psurv( p_int_r( a, t ) );
    assert( psurv > 0.0 );
    const Real target( psurv * rnd );

    const Real thresholdDistance( this->CUTOFF_H * sqrt( 6.0 * D * t ) );

    p_r_params params = { this, t, target };

    gsl_function F;
    if( a <= thresholdDistance )
    {
        F.function = reinterpret_cast<typeof(F.function)>( &p_r_F );
    }
    else
    {
        // p_int_r < p_int_r_free
/*        if( p_int_r_free( a, t ) < target )
        {
            return a;
            }*/

        F.function = reinterpret_cast<typeof(F.function)>( &p_r_free_F );
    }

    F.params = &params;

    const Real low( 0.0 );
    const Real high( a );

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    const Real r( findRoot( F, solver, low, high, 1e-18, 1e-12,
                            "FirstPassageGreensFunction::drawR" ) );
  
    gsl_root_fsolver_free( solver );

    return r;
}



const std::string FirstPassageGreensFunction::dump() const
{
    std::ostringstream ss;
    ss << "D = " << this->getD() << ", a = " << this->geta() << std::endl;
    return ss.str();
}    

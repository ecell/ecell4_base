//#define NDEBUG
//#define BOOST_DISABLE_ASSERTS

#include <iostream>
#include <stdexcept>
#include <vector>
#include <sstream>

#include <boost/bind.hpp>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
//#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sum.h>

#include "factorial.hpp"
#include "bessel.hpp"

//#include "HalfOrderBesselGenerator.hpp"

#include "FirstPassageNoCollisionPairGreensFunction.hpp"



FirstPassageNoCollisionPairGreensFunction::
FirstPassageNoCollisionPairGreensFunction( const Real D ) 
    :
    PairGreensFunction( D, 0, 0 ),
    a( INFINITY )
{
    ; // do nothing
}

FirstPassageNoCollisionPairGreensFunction::
~FirstPassageNoCollisionPairGreensFunction()
{
    ; // do nothing
}

void FirstPassageNoCollisionPairGreensFunction::seta( const Real a )
{
    this->a = a;

//    clearAlphaTable();
}


const Real
FirstPassageNoCollisionPairGreensFunction::p_survival( const Real t,
                                                       const Real r0 ) const
{
    const Real D( getD() );
    const Real a( geta() );

    const Real Dt( D * t );
    const Real a_r( 1.0 / a );

    const Real PIr0( M_PI * r0 );

    const Real angle_factor( PIr0 * a_r );
    const Real exp_factor( - Dt * M_PI * M_PI * a_r * a_r );

    const Real threshold( exp( exp_factor ) * this->TOLERANCE );
    const unsigned int 
        i_max( static_cast<unsigned int>( ceil( a * 
                                                sqrt( log( threshold ) / Dt ) / 
                                                M_PI ) ) );

    Real p( 2.0 * a / PIr0 );
    Real sign( 1.0 );
    unsigned int i( 1 );
    while( true )
    {
        const Real term( sign * 
                         exp( i * i * exp_factor ) * 
                         sin( i * angle_factor ) / i );
        
        p += term;

        if( i >= i_max )
        {
            break;
        }

        sign = -sign;
        ++i;
    }

    return p;
}




const Real
FirstPassageNoCollisionPairGreensFunction::p_int_r( const Real r,
                                                    const Real t,
                                                    const Real r0 ) const
{
    const Real D( getD() );
    const Real a( geta() );

    const Real Dt( D * t );
    const Real a_r( 1.0 / a );

    const Real PIr0( M_PI * r0 );
    const Real PIr( M_PI * r );

    const Real r0_angle_factor( PIr0 * a_r );
    const Real r_angle_factor( PIr * a_r );
    const Real exp_factor( - Dt * M_PI * M_PI * a_r * a_r );

    const Real threshold( exp( exp_factor ) * this->TOLERANCE );
    const unsigned int 
        i_max( static_cast<unsigned int>( ceil( a * 
                                                sqrt( log( threshold ) / Dt ) / 
                                                M_PI ) ) );

    Real p( 2.0  / ( M_PI * PIr0 ) );
    unsigned int i( 1 );
    while( true )
    {
        Real sin_r;
        Real cos_r;
        sincos( r_angle_factor * i, &sin_r, &cos_r );

        const Real term1( exp( exp_factor * i * i ) * 
                          sin( r0_angle_factor * i ) );
        const Real term2( a * sin_r - PIr * i * cos_r );
        const Real term( term1 * term2 / ( i * i ) );
        
        p += term;

        if( i >= i_max )
        {
            break;
        }

        ++i;
    }

    return p;
}




/*
const Real 
FirstPassageNoCollisionPairGreensFunction::drawTime( const Real rnd, 
                                                     const Real r0 ) const
{
   const Real a( this->geta() );

   THROW_UNLESS( std::invalid_argument, rnd <= 1.0 && rnd >= 0.0 );
   THROW_UNLESS( std::invalid_argument, r0 >= 0.0 && r0 <= a );

   if( r0 == a || a == 0.0 )
   {
       return 0.0;
   }

   Real low( 1e-6 );
   Real high( 1.0 );

   p_survival_params params = { this, r0, rnd };

   gsl_function F = 
       {
           reinterpret_cast<typeof(F.function)>( &p_survival_F ),
           &params 
       };

    // adjust high and low to make sure that f( low ) and f( high ) straddle.
    while( GSL_FN_EVAL( &F, high ) < 0.0 )
    {
	high *= 10;
	printf( "drawTime: adjusting high: %g\n", high );
	if( fabs( high ) >= 1e10 )
	{
	    std::cerr << "Couldn't adjust high. F(" << high <<
		") = " << GSL_FN_EVAL( &F, high ) << "; r0 = " << r0 << 
		", " << dump() << std::endl;
	    throw std::exception();
	}
    }

    Real low_value( GSL_FN_EVAL( &F, low ) );
    while( low_value > 0.0 )
    {
	low *= .1;

        const Real low_value_new( GSL_FN_EVAL( &F, low ) );

	printf( "drawTime: adjusting low: %g, F = %g\n", low, low_value_new );

	if( fabs( low ) <= this->MIN_T || 
            fabs( low_value - low_value_new ) < TOLERANCE ) 
	{
	    std::cerr << "Couldn't adjust low.  Returning MIN_T (= "
		      << this->MIN_T << "); F(" << low <<
		") = " << GSL_FN_EVAL( &F, low ) << "; r0 = " << r0 << ", "
		      << dump() << std::endl;
	    return this->MIN_T;
	}

        low_value = low_value_new;
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

	const int status( gsl_root_test_interval( low, high, this->MIN_T, 
						  this->TOLERANCE ) );

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
  
    // printf("%d\n", i );


    Real t( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );

    return t;
}
*/
const Real 
FirstPassageNoCollisionPairGreensFunction::drawR( const Real rnd, 
                                                  const Real r0, 
                                                  const Real t ) const
{
    return 0.0;
}
    
const Real 
FirstPassageNoCollisionPairGreensFunction::drawTheta( const Real rnd,
                                                      const Real r, 
                                                      const Real r0, 
                                                      const Real t ) const
{
    return 0.0;
}


#define HAVE_INLINE

//#define NDEBUG
//#define BOOST_DISABLE_ASSERTS

#include <exception>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_elljac.h>


#include "FirstPassageGreensFunction.hpp"


const Real ellipticTheta4zero( const Real q )
{
  const Int N(10000);
  Real value( 1.0 );

  for( Int n( 1 ); n < N+1; ++n )
    {
      // here I assume this series converges rapidly and
      // n doesn't become larger than, say, 10.
      const Real qtothe2nm1( gsl_pow_int( q, 2 * n - 1 ) );

      Real term1( 1.0 - qtothe2nm1 * q );
      Real term2( 1.0 - qtothe2nm1 );
      term2 *= term2;

      const Real term( term1 * term2 );
      const Real value_prev( value );
      value *= term;
      
      if( fabs( value - value_prev ) < 1e-8 )
	{
	  break;
	}

      value *= term;
    }

  return value;
}

const Real 
FirstPassageGreensFunction::p_survival( const Real time ) const
{
  const Real D( getD() );
  const Real a( geta() );
  const Real asq( a * a );
  const Real PIsq( M_PI * M_PI );

  const Real q( - ( D * PIsq * time ) / asq );

  return 1.0 - ellipticTheta4zero( exp( q ) );
} 

const Real 
FirstPassageGreensFunction::drawExitTime( const Real rnd, const Real r,
					  const Real time ) const
{

}


const Real 
FirstPassageGreensFunction::drawR( const Real rnd, const Real r,
				   const Real t ) const
{

}

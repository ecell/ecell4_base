#define HAVE_INLINE

//#define NDEBUG
//#define BOOST_DISABLE_ASSERTS

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


#include "FirstPassageGreensFunction.hpp"



//FIXME: move to better place
const Real CUTOFF( 1e-10 );


/*
  EllipticTheta[4,0,q]
*/
static const Real ellipticTheta4zero( const Real q )
{
  assert( q < 1.0 );

  const Int N( 100 );
  Real value( 1.0 );

  Real q_n( q );
  Real q_2n( 1.0 );

  for( Int n( 1 ); n <= N; ++n )
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
	  //printf("%d\n",n);
	  return value;
	}

      q_n *= q;  // q_(++n)
    }

  std::cerr << "WARNING: ellipticTheta4zero: didn't converge." << std::endl;
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

  return 1.0 - ellipticTheta4zero( exp( q ) );
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
FirstPassageGreensFunction::p_r_int( const Real r, const Real t ) const
{
  Real value( 0.0 );

  const Real p_free( this->p_free_int( r, t ) );

  // p_r_int is always smaller than p_free.
  if( fabs( p_free ) < CUTOFF )
    {
      return 0.0;
    }

  const Real D( getD() );
  const Real a( geta() );
  const Real asq( a * a );
  const Real PIsq( M_PI * M_PI );

  const Real PIr( M_PI * r );
  const Real PIr_a( PIr / a );
  const Real DtPIsq_asq( D * t * PIsq / asq );

  const Real factor( 2.0 / ( a * M_PI ) );
  const Real p_free_factor( p_free * factor );

  const Int N( 1000 );
  long int n( 1 );
  while( true )
    {
      const Real term1( exp( - n * n * DtPIsq_asq ) );
      
      const Real term2( a * sin( n * PIr_a ) );
      const Real term3( n * PIr * cos( n * PIr_a ) );

      const Real term( term1 * ( term2 - term3 ) / n );
      value += term;

      printf("%ld %g %g %g %g %g %g\n", n, p_free_factor,
	     value, term, term1, term2, term3 );

      if( fabs( value ) * 1e-10 >= fabs( term ) )
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

  printf( "Dt/a^2 %g\tfree %g\n", D*t/(a*a),p_free_int( r, t ) );
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

  const Int N( 100 );

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
	  std::cerr << "p_r_int: didn't converge; " << n << " " << value 
		    << std::endl;
	  break;
	}

      ++n;
    }

  const Real factor( 1.0 / ( sqrt( 2 ) * PIsq * pow(D * t, 1.5) ) );

  return value * factor;
} 




const Real 
FirstPassageGreensFunction::drawExitTime( const Real rnd, const Real r,
					  const Real t ) const
{

}


const Real 
FirstPassageGreensFunction::drawR( const Real rnd, const Real r,
				   const Real t ) const
{

}



#ifdef __FPGF_TEST__

int main()
{
  printf("%g\n",ellipticTheta4zero( 0.9 ) );
}



#endif

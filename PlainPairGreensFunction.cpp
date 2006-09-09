#define HAVE_INLINE

//#define NDEBUG
#define BOOST_DISABLE_ASSERTS

#include <exception>
#include <vector>

#include <boost/array.hpp>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_roots.h>


#include "HalfOrderBesselGenerator.hpp"

#include "bessel.c"

#include "PlainPairGreensFunction.hpp"



PlainPairGreensFunction::PlainPairGreensFunction( const Real D, 
						  const Real kf, 
						  const Real Sigma )
  :
  PairGreensFunction( D, kf, Sigma ),
  kD( 4.0 * M_PI * getSigma() * getD() ),
  alpha( ( 1.0 + ( getkf() / getkD() ) ) * ( sqrt( getD() ) / getSigma() ) )
{
  ; // do nothing
}

PlainPairGreensFunction::~PlainPairGreensFunction()
{
  ; // do nothing
}





const Real 
PlainPairGreensFunction::p_corr_R2( const Real u, 
				    const p_corr_R2_params* const params )
{
  const Real SIGMA2KFp1( 2.0 * params->Sigma * params->kf + 1.0 );
  const Real SIGMA2U( 2.0 * u * params->Sigma );

  const Real costheta( cos( params->theta ) );

  const int order_min( params->order );
  const int order_max( params->order_max );
  HalfOrderBesselGenerator us( u * params->Sigma, order_min, order_max );
  HalfOrderBesselGenerator ur( u * params->r, order_min, order_max );
  HalfOrderBesselGenerator ur0( u * params->r0, order_min, order_max );

  Real result( 0.0 );
  for( int order( order_min ); order <= order_max; ++order )
    {
      const Real js( us.j( order ) );
      const Real ys( us.y( order ) );
      const Real jps( us.jp( order ) );
      const Real yps( us.yp( order ) );
      const Real jr( ur.j( order ) );
      const Real yr( ur.y( order ) );
      const Real jr0( ur0.j( order ) );
      const Real yr0( ur0.y( order ) );

    
      // R1 / ( R1^2 + R2^2 ) * ( R1 F2 + R2 F2 )
      
      // below I rewrote the equation in the way that 
      // (1) j and y are multiplied first, before j^2 or y^2 occurs,
      //     of which absolute values of exponents can be huge 
      //    (either positive or negative).
      // (2) avoid roundoff error caused by j? - j?.
      //     (including j? - jp? etc.)  assumed y? - y? is ok (?).
      
      const Real R1a( SIGMA2KFp1 * js );
      const Real R1b( SIGMA2U * jps );
      
      const Real R1F1( ( ( R1a * jr ) * jr0 - ( R1a * yr ) * yr0 ) + 
		       ( R1b * yr ) * yr0 - ( R1b * jr ) * jr0 );

      const Real R2( SIGMA2KFp1 * ys - SIGMA2U * yps ); 
      const Real F2( jr * yr0 + jr0 * yr );

      const Real R1F1_plus_R2F2( R1F1 + R2 * F2 );
      const Real num( R1a * R1F1_plus_R2F2 - R1b * R1F1_plus_R2F2 );

      // R1^2 + R2^2, roundoff error here wouldn't be a big problem, though.
      const Real den( R2 * R2 + R1a * R1a - 2.0 * R1a * R1b + R1b * R1b );

      const Real lgnd( gsl_sf_legendre_Pl( order, costheta ) );
      const Real coeff_corr( order + order + 1.0 );

      /*
	const double _R1( SIGMA2KFp1 * js - SIGMA2U * jps );
	const double _R2( SIGMA2KFp1 * ys - SIGMA2U * yps );
  
	const double _R1F1( jr*jr0*_R1 - (_R1*yr)*yr0 );
	const double _F2( yr0*jr + jr0*yr );
      
	const double result2( (_R1 * _R1F1+ (_R1*_R2)*_F2) / (_R1*_R1+_R2*_R2) );
      */

      result += (num/den) * lgnd * coeff_corr;

      //      printf("%g %g\n",((num/den)-result2),result2);

    }

  const Real exp_term( exp( - params->D * u * u * params->t ) * u );
  result *= exp_term;




  if( isnan( result ) )
    {
      std::cerr << "NaN in p_corr_R" << std::endl;
      std::cerr << u << ' ' << order_min << ' ' << exp_term << std::endl;
      //      std::cerr << "R1F1 " << R1F1 << " R2 F2" << R2 << ' ' << F2 << " result" << result <<std::endl;
      throw std::exception(); //("NaN in p_corr_R");
    }

  return result;
}

const Real 
PlainPairGreensFunction::p_corr_R( const Real u, 
				   const p_corr_R_params* const params )
{
  Real result;

  const Real SIGMA2KFp1( params->Sigma * 2 * params->kf + 1.0 );
  const Real SIGMA2U( params->Sigma * 2 * u );

  const int order( params->order );

  Real jr,yr,jr0,yr0,js,ys,jps,yps,tmp_;
  bessjy(u*params->Sigma,order+0.5,&js,&ys,&jps,&yps);
  bessjy(u*params->r,order+0.5,&jr,&yr,&tmp_,&tmp_);
  bessjy(u*params->r0,order+0.5,&jr0,&yr0,&tmp_,&tmp_);

  // R1 / ( R1^2 + R2^2 ) * ( R1 F2 + R2 F2 )
  
  // below I rewrote the equation in the way that 
  // (1) j and y are multiplied first, before j^2 or y^2 occurs,
  //     of which absolute values of exponents can be huge 
  //    (either positive or negative).
  // (2) avoid roundoff error caused by j? - j?.
  //     (including j? - jp? etc.)  assumed y? - y? is ok (?).
  
  const Real R1a( SIGMA2KFp1 * js );
  const Real R1b( SIGMA2U * jps );
  
  const Real R1F1( ( ( R1a * jr ) * jr0 - ( R1a * yr ) * yr0 ) + 
		   ( R1b * yr ) * yr0 - ( R1b * jr ) * jr0 );
  
  const Real R2( SIGMA2KFp1 * ys - SIGMA2U * yps ); 
  const Real F2( jr * yr0 + jr0 * yr );
  
  const Real R1F1_plus_R2F2( R1F1 + R2 * F2 );
  const Real num( R1a * R1F1_plus_R2F2 - R1b * R1F1_plus_R2F2 );
  
  // R1^2 + R2^2, roundoff error here wouldn't be a big problem, though.
  const Real den( R2 * R2 + R1a * R1a - 2.0 * R1a * R1b + R1b * R1b );

  /*
    const double _R1( SIGMA2KFp1 * js - SIGMA2U * jps );
    const double _R2( SIGMA2KFp1 * ys - SIGMA2U * yps );
  
    const double _R1F1( jr*jr0*_R1 - (_R1*yr)*yr0 );
    const double _F2( yr0*jr + jr0*yr );
  
    const double num( _R1 * _R1F1+ (_R1*_R2)*_F2 );
    const double den( _R1*_R1+_R2*_R2 );
  */
  
  const Real exp_term( exp( - params->D * u * u * params->t ) * u );

  result = ( num / den ) * exp_term;

  if( isnan( result ) )
    {
      std::cerr << "Error: NaN in p_corr_R" << std::endl;
    }

  return result;
}



const Real 
PlainPairGreensFunction::p_corr( const Real r, const Real r0, 
				 const Real theta, const Real t ) const
{
  //  int NU=(int)(100*pow((D/0.05*time),-0.3)); 
  //Real umax=25/sqrt(D/0.05*time); 
    

  // the point that the scaling factor becomes exp( -30 ) ~= 1e-13
  const Real umax( sqrt( 30.0 / ( getD() * t ) ) ); 
    
  const int maxorder( 100 );
    
  // initially calculate orders \in [0,3], use of recursion implemented in
  // HalfBesselGenerator is justified at least two orders, plus overhead.
  int order_step( 4 );

  gsl_function p_corr_R2_F;
  p_corr_R2_F.function = 
    reinterpret_cast<typeof(p_corr_R2_F.function)>( &p_corr_R2 );

  //              params = { nmin, nmax-1, r, r0, theta, t, Sigma, D, kf };
  p_corr_R2_params params = { 0, order_step-1, r, r0, theta, t, getSigma(), 
			      getD(), getkf() };
  p_corr_R2_F.params = reinterpret_cast<typeof(p_corr_R2_F.params)>( &params );

  gsl_integration_workspace* 
    workspace( gsl_integration_workspace_alloc( 1000 ) );

  Real inttot( 0.0 );
  while(true)
    {			
      Real integral;
      Real error;
	
      // abs_err >> 1 because the integral can be huge.
      // instead use rel_err.
      gsl_integration_qag( &p_corr_R2_F, 0.0, 
			   umax, 
			   1e4, 1e-18, 1000,
			   GSL_INTEG_GAUSS61,
			   workspace, &integral, &error );
      
      //      printf("%g %g\n",integral,error);
      inttot += integral;

      //printf("%d %g %g\n",params.order,inttot, integral );
	
      // truncate when converged enough.
      const Real integral_abs( fabs( integral ) );
      if( integral_abs / fabs( inttot ) < 1e-10  )
	{
	  break;
	}
	
      if( params.order >= maxorder )
	{
	  std::cerr << "Didn't converge. " << params.order << std::endl;
	  // throw std::exception();//"Didn't converge.");
	  break;
	}
	

      params.order = params.order_max + 1;
      params.order_max = params.order + order_step;

      if( order_step <= 10 )
	{
	  ++order_step;
	}
      //printf("%d %d\n", params.order,params.order_max);
    }	

  gsl_integration_workspace_free( workspace );

  //printf("%d\n",params.nmax);
    
  inttot /= - ( M_PI*4.0*sqrt( r * r0 ) );

  //printf("%d\n",params.order_max);
    
  return inttot;
}
  

const Real 
PlainPairGreensFunction::p_free( const Real r, const Real r0, 
				 const Real theta, const Real t ) const
{
  const Real Dt4( 4.0 * getD() * t );
  const Real Dt4PI( Dt4 * M_PI );

  return exp( ( 2.0 * r * r0 * cos( theta ) - r * r - r0 * r0 ) / Dt4 )
    / sqrt( Dt4PI * Dt4PI * Dt4PI );
}

const Real 
PlainPairGreensFunction::p_tot( const Real r, const Real r0, 
				const Real theta, const Real t ) const
{
  const Real factor( 2.0 * M_PI * r * r * sin( theta ) );

  const Real p_free( this->p_free( r, r0, theta, t ) * factor );

  if( p_free <= PlainPairGreensFunction::P_CUTOFF )
    {
      return 0.0;
    }

  const Real p_corr( this->p_corr( r, r0, theta, t ) * factor );

  Real p_tot( p_free + p_corr );
		    
						
  if( p_tot < 0.0 )
    {
      if( - p_tot < P_CUTOFF )
	{
	  std::cerr << "p_tot is negative, but" <<
	    " abs( p_tot ) < EPSILON." <<
	    "setting zero.\n"
		    << "p_tot = " << p_tot
		    << ", p_corr = " << p_corr
		    << ", p_free = " << p_free
		    << std::endl;
	  std::cerr << "t " << t << "\tr0 " << r0 << "\tr " << r
		    << "\ttheta " << theta << std::endl;
	  //	    exit(1);
	  p_tot = 0.0;
	}
      else
	{
	  std::cerr << "WARNING: p_tot is negative, and" <<
	    " abs( p_tot ) >= EPSILON. Setting zero.\n"
		    << "p_tot = " << p_tot
		    << ", p_corr = " << p_corr
		    << ", p_free = " << p_free
		    << std::endl;
	  std::cerr << "t " << t << "\tr0 " << r0 << "\tr " << r
		    << "\ttheta " << theta << std::endl;

	  p_tot = 0.0;
	  //exit( 1 );
	}
    }
    
  return p_tot;
}


//  virtual const Real p_t_given_r0( const Real t, const Real r0 );

/**
   Calculates exp( x^2 ) * erfc( x )

   See asymptotic expansion here:
   http://en.wikipedia.org/wiki/Error_function
*/  
inline const Real expxsq_erfc( const Real x )
{
  Real result;

  const Real xsq( x * x );
  if( x > 26.0 )
    {
      static const Real M_1_SQRTPI( M_2_SQRTPI * 0.5 ); 

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
inline const Real W( const Real a, const Real b )
{

  // exp( 2 a b + b^2 ) erfc( a + b ) == 
  //               exp( - a^2 ) exp( ( a + b )^2 ) erfc( a + b )
  return exp( - a * a ) * expxsq_erfc( a + b );
}

const Real 
PlainPairGreensFunction::p_irr_radial( const Real r,
				       const Real t,
				       const Real r0 ) const
{
  //  printf("irrp %g %g %g\n",r,r0,t);
  const Real sqrtD( sqrt( getD() ) );
  const Real alpha( this->getalpha() );

  const Real Dt4( 4.0 * this->getD() * t );

  const Real r_plus_r0_minus_2sigma( r + r0 - 2.0 * this->getSigma() );

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
PlainPairGreensFunction::p_reaction( const Real tsqrt, const Real r0 ) const
{
  const Real kD( this->getkD() );
  const Real kf( this->getkf() );
  const Real Sigma( this->getSigma() );
  const Real D( this->getD() );
  const Real alpha( this->getalpha() );

  const Real sqrtD( sqrt( D ) );

  const Real r0_m_Sigma_over_sqrt4D_t( ( r0 - Sigma ) 
				       / ( ( sqrtD + sqrtD ) * tsqrt ) );

  const Real Wf( W( r0_m_Sigma_over_sqrt4D_t, alpha * tsqrt ) );
  const Real factor( Sigma * kf / ( r0 * ( kf + kD ) ) );

  return factor * ( erfc( r0_m_Sigma_over_sqrt4D_t ) - Wf );
}

const Real 
PlainPairGreensFunction::p_reaction_deriv( const Real tsqrt, 
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
PlainPairGreensFunction::p_reaction_fdf( const Real tsqrt, 
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




const Real 
PlainPairGreensFunction::p_reaction_F( const Real tsqrt, 
				       const p_reaction_params* const params )
{
  const PlainPairGreensFunction* const gf( params->gf ); 
  const Real r0( params->r0 );
  const Real rnd( params->rnd );

  return gf->p_reaction( tsqrt, r0 ) - rnd;
}

const Real 
PlainPairGreensFunction::
p_reaction_deriv_F( const Real tsqrt, const p_reaction_params* const params )
{
  const PlainPairGreensFunction* const gf( params->gf ); 
  const Real r0( params->r0 );
  //  const Real rnd( params->rnd );

  return gf->p_reaction_deriv( tsqrt, r0 );
}

void
PlainPairGreensFunction::
p_reaction_fdf_F( const Real tsqrt, const p_reaction_params* const params,
		  Real* const f, Real* const df )
{
  const PlainPairGreensFunction* const gf( params->gf ); 
  const Real r0( params->r0 );
  const Real rnd( params->rnd );

  gf->p_reaction_fdf( tsqrt, r0, f, df );
  *f -= rnd;
}


const Real PlainPairGreensFunction::drawTime( const Real rnd, 
					      const Real r0, 
					      const Real maxt ) const
{
  assert( rnd <= 1.0 && rnd >= 0.0 );

  {
    const Real maxp( p_reaction( maxt, r0 ) );

    if( rnd >= maxp )
      {
	return INFINITY;
      }
  }

  p_reaction_params params = { this, r0, rnd };
  gsl_function_fdf FDF = 
    {
      reinterpret_cast<typeof(FDF.f)>( &p_reaction_F ),
      reinterpret_cast<typeof(FDF.df)>( &p_reaction_deriv_F ),
      reinterpret_cast<typeof(FDF.fdf)>( &p_reaction_fdf_F ),
      &params 
    };

  // This initial guess must be smaller than the answer, or
  // newton-type iteration can undershoot into tsqrt < 0, where
  // value of p_reaction is undefined.
  const Real initialGuess( 1e-100 );

  //const gsl_root_fdfsolver_type* solverType( gsl_root_fdfsolver_steffenson );
  const gsl_root_fdfsolver_type* solverType( gsl_root_fdfsolver_secant );
  gsl_root_fdfsolver* solver( gsl_root_fdfsolver_alloc( solverType ) );
  gsl_root_fdfsolver_set( solver, &FDF, initialGuess );

  const Index maxIter( 100 );

  Real tsqrt( initialGuess );
  Index i( 0 );
  while( true )
    {
      gsl_root_fdfsolver_iterate( solver );
      const Real tsqrt_prev( tsqrt );
      tsqrt = gsl_root_fdfsolver_root( solver );

      int status( gsl_root_test_delta( tsqrt, tsqrt_prev, 1e-18, 1e-12 ) );

      if( status == GSL_CONTINUE )
	{
	  if( i >= maxIter )
	    {
	      gsl_root_fdfsolver_free( solver );
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
  
  gsl_root_fdfsolver_free( solver );

  return gsl_pow_2( tsqrt );
} 




const Real PlainPairGreensFunction::drawR( const Real rnd, 
					   const Real r0, 
					   const Real t ) const
{
  Real r;

  const Index tableSize( 100 );

  boost::array<Real, tableSize> pTable;
  pTable[0] = 0.0;

  const Real minR( this->getMinR() );
  const Real maxR( this->getMaxR( t ) );
  const Real rStep( ( maxR - minR ) / ( tableSize - 1) );

  Real r_prev( p_irr_radial( minR, t, r0 ) );
  Index i( 1 );
  for( ; i < tableSize; ++i )
    {
      const Real r( minR + i * rStep );
      const Real p( p_irr_radial( r, t, r0 ) );
      const Real value( ( r_prev + p ) * 0.5 );
      pTable[i] = pTable[i-1] + value;

      if( value < pTable[i] * std::numeric_limits<double>::epsilon() ) 
	{
	  break;   // truncation; pTable is valid in [0,i].
	}

      r_prev = p;
    }

  if( i >= tableSize-1 )
    {
      std::cerr << "drawR: p didn't converge within the range of table." <<
	std::endl;
    }

  const Real targetPoint( rnd * pTable[i] );
  const size_t lowerBound( gsl_interp_bsearch( pTable.data(), targetPoint, 
					       0, i ) );

  const Real low( minR + lowerBound * rStep );
      
  if( pTable[lowerBound+1] - pTable[lowerBound] != 0.0 )
    {
      r = low + rStep * ( targetPoint - pTable[lowerBound] ) / 
	( pTable[lowerBound+1] - pTable[lowerBound] );
    }
  else
    {
      // this can happen when rnd is or is too close to 1.0.
      r = low;
    }

  return r;
}


const Real PlainPairGreensFunction::
Rn( const Int order, const Real r, const Real r0, const Real t,
    gsl_integration_workspace* workspace ) const
{
  Real integral;
  Real error;

  gsl_function p_corr_R_F;
  p_corr_R_F.function = 
    reinterpret_cast<typeof(p_corr_R_F.function)>( &p_corr_R );

  p_corr_R_params params = { order, r, r0, t, getSigma(), getD(), getkf() };
  p_corr_R_F.params = reinterpret_cast<typeof(p_corr_R_F.params)>( &params );


  const Real umax( sqrt( 30.0 / ( this->getD() * t ) ) ); 
  // abs_err >> 1 because the integral can be huge.
  // instead use rel_err.
  gsl_integration_qag( &p_corr_R_F, 0.0, 
		       umax,
		       1e4, 1e-18, 1000,
		       GSL_INTEG_GAUSS61,
		       workspace, &integral, &error );
  
  return integral;
}


const Real PlainPairGreensFunction::
p_corr_RnTable( const Real theta, const Real r, const Real r0,
		const Real t, const RealVector& RnTable )
{
  Real result( 0.0 );

  const Index tableSize( RnTable.size() );
  
  RealVector PnTable( tableSize );
  gsl_sf_legendre_Pl_array( tableSize-1, cos( theta ), PnTable.data() );

  for( Index order( 0 ); order < tableSize; ++order )
    {
      const Real Cn( 2.0 * order + 1 );

      //      printf( "p %d %g\n", order, PnTable[order] );

      result += Cn * PnTable[order] * RnTable[order];
    }

  result /= 4.0 * M_PI * sqrt( r * r0 );

  return result;
}

const Real PlainPairGreensFunction::
p_tot_RnTable( const Real theta, const Real r, const Real r0,
	       const Real t, const RealVector& RnTable ) const
{
  const Real factor( 2.0 * M_PI * r * r * sin( theta ) );

  const Real p_corr( this->p_corr_RnTable( theta, r, r0, t, RnTable ) ); 
  const Real p_free( this->p_free( r, r0, theta, t ) );

  return ( p_free + p_corr ) * factor;
}

const Real PlainPairGreensFunction::drawTheta( const Real rnd,
					       const Real r, 
					       const Real r0, 
					       const Real t ) const
{
  Real theta;

  RealVector RnTable;
  RnTable.reserve( 8 );

  const Int MAXORDER( 100 );

  gsl_integration_workspace* 
    workspace( gsl_integration_workspace_alloc( 1000 ) );

  Real RnCum( 0.0 );
  for( Int order( 0 ); order <= MAXORDER; ++order )
    {
      const Real Rn( this->Rn( order, r, r0, t, workspace ) );
      RnTable.push_back( Rn );
      RnCum += Rn;
      // truncate when converged enough.
      if( fabs( Rn ) < 
	  fabs( RnCum * std::numeric_limits<double>::epsilon() ) )
	{
	  printf("%d\n",order );
	  break;
	}
    }

  gsl_integration_workspace_free( workspace );

  const Index tableSize( 200 );
  const Real thetaStep( M_PI / (tableSize-1) );

  boost::array<Real, tableSize> pTable;
  pTable[0] = 0.0;

  Real p_prev( this->p_tot_RnTable( 0.0, r, r0, t, RnTable ) );
  Index i( 1 );
  for( ; i < tableSize; ++i )
    {
      const Real theta( thetaStep * i );

      Real p( this->p_tot_RnTable( theta, r, r0, t, RnTable ) );
      if( p < 0.0 )
	{
	  p = 0.0;
	}

      const Real value( ( p_prev + p ) * 0.5 );

      pTable[i] = pTable[i-1] + value;

      if( value < pTable[i] * std::numeric_limits<Real>::epsilon() ) 
	{
	  break;   // truncation; pTable is valid in [0,i].
	}

      p_prev = p;
    }

  // Ideally, p_irr_r below can work as a normalization factor,
  const Real p_irr_r( this->p_irr_radial( r, t, r0 ) );
  const Real relerror( (pTable[i]*thetaStep- p_irr_r)/p_irr_r );
  //  printf("%g %g\n", pTable[i-1]*thetaStep, p_irr_r );
  if( fabs( relerror ) >= 1e-2 )
    {
      std::cerr << "drawTheta: relative error estimate is large: " << relerror 
		<< "( r= " << r << ", r0= " << r0 << ", t= " << t 
		<< ")." << std::endl;
    }

  if( i >= tableSize-1 )
    {
      std::cerr << "drawTheta: p didn't converge within the range of table." <<
	std::endl;
    }

  const Real targetPoint( rnd * pTable[i] );
  const size_t lowerBound( gsl_interp_bsearch( pTable.data(), targetPoint, 
					       0, i ) );
  const Real low( lowerBound * thetaStep );
      
  if( pTable[lowerBound+1] - pTable[lowerBound] != 0.0 )
    {
      theta = low + thetaStep * ( targetPoint - pTable[lowerBound] ) / 
	( pTable[lowerBound+1] - pTable[lowerBound] );
    }
  else
    {
      // this can happen when rnd is or is too close to 1.0.
      theta = low;
    }


  return theta;
}









#ifdef __TEST_PLAINPAIRGREENSFUNCTION

int main(int n, char *input[])
{
  const Real sigma( atof(input[1]) );
  const Real D( atof(input[2]) );
  const Real kf( atof(input[3]) );

  static const Real H = 4.0;
  static const int NR_FS = 20;
  static const int NR0_FS = 20;
  static const int NTIME_FS = 20;
  static const int NTHETA_FS = 39;
	

  PlainPairGreensFunction gf( D, kf, sigma );

  //double theta_fs[NR_FS][NTHETA_FS][NR0_FS][NTIME_FS],time_fs[NTIME_FS],r0_fs[NR0_FS][NTIME_FS],r_fs[NR_FS][NR0_FS][NTIME_FS];
	
  Real4DArray p_tot( boost::extents[NR_FS][NTHETA_FS][NR0_FS][NTIME_FS] );
  Real4DArray P_tot( boost::extents[NR_FS][NTHETA_FS][NR0_FS][NTIME_FS] );
  Real3DArray p_tot_r( boost::extents[NR_FS][NR0_FS][NTIME_FS] );
  Real4DArray p_tot_theta( boost::extents
			   [NR_FS][NTHETA_FS][NR0_FS][NTIME_FS] );


  
  for( int i( 0 ); i <= 1000000 ; ++i )
    {
      Real result= gf.drawR( 0.5, 1e-8, 1e-6 );
    }

  exit(1);

  //  double a = gf.p_tot( 1.32072e-8, 1.998e-8, 0.668163, 1e-6);
  //  exit(1);
  //  Set_Parameter_Range( sigma, D, kf );

  //t 0.001        r0 3.09939e-07  r 3.9502e-07    theta 1.24678

  for(int i=0;i<NR_FS;i++)
    {
      for(int j=0;j<NTHETA_FS;j++)
	{
	  for(int k=0;k<NR0_FS;k++)		
	    {
	      std::cerr << k << std::endl;
	      for(int l=0;l<NTIME_FS;l++)
		{
		  const Real t( gf.getTByScale( Real( l ) / NTIME_FS ) );
		  const Real r0( gf.getR0ByScale( Real( k ) / NR0_FS, t ) );
		  const Real r( gf.getRByScale( Real( i ) / NR_FS, t ) );
		  const Real theta( gf.getThetaByScale( Real( j ) / NTHETA_FS,
							t, r0, r ) );


		  //std::cerr << "t " << t << "\tr0 " << r0 << "\tr " << r
		  //			    << "\ttheta " << theta << std::endl;
		  p_tot[i][j][k][l] = gf.p_tot( r, r0, theta, t );
		  //		  std::cerr << p_tot[i][j][k][l] << std::endl;  

		}
	    }
	}	
    }		

#if 0
  //	BUILDING p(r|r0,t), INTEGRATING OUT THETA
  for(int i=0;i<NR_FS;i++)
    {
      for(int j=1;j<NTHETA_FS;j++)
	{
	  for(int k=0;k<NR0_FS;k++)
	    {
	      for(int l=0;l<NTIME_FS;l++)
		{
		  p_tot_r[i][k][l] +=
		    ( p_tot[i][j-1][k][l] + p_tot[i][j][k][l] ) / 2 *
		    ( theta_fs[i][j][k][l] - theta_fs[i][j-1][k][l] );					
		  //					printf("%.12f\t%f\t%f\t%f\n",p_tot_r[i][k][l],time_fs[l],r0_fs[k][l],r_fs[i][k][l]);
		}
	    }
	}
    }

  // BUILDING p(theta|r,r0,t)=p(r,theta|r0,t)/p(r|r0,t)							
  for(int i=0;i<NR_FS;i++)
    {
      for(int j=0;j<NTHETA_FS;j++)
	{
	  for(int k=0;k<NR0_FS;k++)
	    {
	      for(int l=0;l<NTIME_FS;l++)
		{
		  p_tot_theta[i][j][k][l] =
		    p_tot[i][j][k][l] / p_tot_r[i][k][l];
		  //					printf("%.16f\t%.16f\t%.16f\n",p_tot[i][j][k][l],p_tot_r[i][k][l],p_tot_theta[i][j][k][l]);
		}
	    }
	}
    }


  char filename1[100],filename2[100],filename3[100],
    filename4[100],filename5[100];
  FILE *prob,*err,*cum,*param,*pfree;

  sprintf(filename1,"../tables/FSProb_D_%.0f_kf_%.0f_tmax_50.dat",D,kf);
  sprintf(filename2,"../tables/FSPfree_D_%.0f_kf_%.0f_tmax_50.dat",D,kf);
  sprintf(filename3,"../tables/Errors_D_%.0f_kf_%.0f_tmax_50.dat",D,kf);
  sprintf(filename4,"../tables/FSCum_D_%.0f_kf_%.0f_tmax_50.dat",D,kf);
  sprintf(filename5,"../tables/FSParam_D_%.0f_kf_%.0f_tmax_50.dat",D,kf);
		
  //	prob=fopen(filename1,"w");
  //	pfree=fopen(filename2,"w");
  err=fopen(filename3,"w");
  cum=fopen(filename4,"w");
  param=fopen(filename5,"w");
	


	
  // BUILDING P(theta|r,r0,t), CUMULATIVE OF p(theta|r,r0,t)
  for(int i=0;i<NR_FS;i++)
    {
      for(int k=0;k<NR0_FS;k++)		
	{
	  for(int l=0;l<NTIME_FS;l++)
	    {
	      P_tot[i][0][k][l]=0;			
	      //				printf("%.12f\t%f\t%f\n",P_tot[i][0][k][l],theta_fs[i][0][k][l],p_tot[i][0][k][l]);
	      fprintf(cum,"%.12f\n",P_tot[i][0][k][l]);
	      fprintf(param,"%f\t%f\t%f\t%f\n",
		      time_fs[l],r0_fs[k][l],r_fs[i][k][l],
		      theta_fs[i][0][k][l]);												
				
	      for(int j=1;j<NTHETA_FS;j++)
		{
		  P_tot[i][j][k][l] = P_tot[i][j-1][k][l] + 
		    ( p_tot_theta[i][j-1][k][l] + p_tot_theta[i][j][k][l] ) / 2
		    * ( theta_fs[i][j][k][l] - theta_fs[i][j-1][k][l] );
					
		  if( isnan( P_tot[i][j][k][l] ) )
		    {
		      P_tot[i][j][k][l]=1.0;
		    }

		  //					printf("%.12f\t%f\t%f\n",P_tot[i][j][k][l],theta_fs[i][j][k][l],p_tot[i][j][k][l]);
		  fprintf(cum,"%.12f\n",P_tot[i][j][k][l]);
		  fprintf(param,"%f\t%f\t%f\t%f\n",
			  time_fs[l],r0_fs[k][l],r_fs[i][k][l],
			  theta_fs[i][j][k][l]);
		}
	      //					printf("\n");
	    }
	}
    }

  //	fclose(prob);	fclose(pfree);
  fclose(err); fclose(cum); fclose(param); 
	
#endif	
  return 0;
}



#endif // __TEST_PLAINPAIRGREENSFUNCTION

//#define NDEBUG
//#define BOOST_DISABLE_ASSERTS

#include <exception>
#include <vector>

#include <boost/array.hpp>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
//#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_roots.h>

#include "bessel.hpp"

#include "freeFunctions.hpp"

#include "HalfOrderBesselGenerator.hpp"

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
    
    const Real SIGMA2KFp1( params->Sigma * 2.0 * params->kf + 1.0 );
    const Real SIGMA2U( params->Sigma * 2.0 * u );
    
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

    assert( ! isnan( result ) );

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
	gsl_integration_qag( &p_corr_R2_F, 0.0, umax, 
			     10, // FIXME: this needs to be adaptively given.
			     1e-8, 
			     1000,
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


const Real 
PlainPairGreensFunction::p_reaction( const Real t, const Real r0 ) const
{
    const Real sqrtt( sqrt( t ) );
    const Real kf( getkf() );
    const Real D( getD() );
    const Real sigma( getSigma() );
    const Real alpha( getalpha() );
    const Real kD( getkD() );


    return __p_reaction_irr( sqrtt, r0, kf, D, sigma, alpha, kD  );
}


const Real 
PlainPairGreensFunction::p_reaction_F( const Real tsqrt, 
				       const p_reaction_params* const params )
{
    const PlainPairGreensFunction* const gf( params->gf ); 
    const Real kf( gf->getkf() );
    const Real D( gf->getD() );
    const Real sigma( gf->getSigma() );
    const Real alpha( gf->getalpha() );
    const Real kD( gf->getkD() );

    const Real r0( params->r0 );
    const Real rnd( params->rnd );

    return __p_reaction_irr( tsqrt, r0, kf, D, sigma, alpha, kD  ) - rnd;
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

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &p_reaction_F ),
	    &params 
	};

    Real low( 1e-100 );
    Real high( sqrt( maxt ) );

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

    return gsl_pow_2( low );
} 




const Real PlainPairGreensFunction::drawR( const Real rnd, 
					   const Real r0, 
					   const Real t ) const
{
    Real r;

    const Index tableSize( 100 );

    boost::array<Real, tableSize> pTable;
    pTable[0] = 0.0;

    // Range of r in this function = 
    //     [ max( Sigma, r0 - H sqrt( 6 D t ), r0(t) + H sqrt( 6 D t ) ];

    const Real Hsqrt6Dt( ( this->H + 1 ) * sqrt( 6.0 * getD() * t ) );

    const Real minR( std::max( getSigma(), r0 - Hsqrt6Dt ) );
    const Real maxR( ( r0 + Hsqrt6Dt ) * 1.1 ); // need to fine tune these

    const Real rStep( ( maxR - minR ) / ( tableSize - 1) );

    Real r_prev( __p_irr( minR, t, r0, getkf(),
                          getD(), getSigma(), getalpha() ) );
    Index i( 1 );

    while( true )
    {
	const Real r( minR + i * rStep );
	const Real p( __p_irr( r, t, r0, getkf(),
                               getD(), getSigma(), getalpha() ) );
	const Real value( ( r_prev + p ) * 0.5 );
	pTable[i] = pTable[i-1] + value;


	// this criterion needs to be decided based upon H.
	if( value < pTable[i] * std::numeric_limits<double>::epsilon() )
	{
	    break;
	}

	if( i >= tableSize - 1 ) 
	{
	    std::cerr << "drawR: p didn't converge within the range of table." <<
		std::endl;
	    break;    
	}

	r_prev = p;

	++i;
    }

    // pTable is valid in [0,i].

    const Real targetPoint( rnd * pTable[i] );
    const size_t lowerBound( gsl_interp_bsearch( &pTable[0], targetPoint, 
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


const Real 
PlainPairGreensFunction::Rn( const Integer order, const Real r, const Real r0,
			     const Real t,
			     gsl_integration_workspace* const workspace,
			     const Real err ) const
{
    Real integral;
    Real error;

    p_corr_R_params params = { order, r, r0, t, getSigma(), getD(), getkf() };
    gsl_function p_corr_R_F = 
	{
	    reinterpret_cast<typeof(p_corr_R_F.function)>( &p_corr_R ),
	    &params
	};

    const Real umax( sqrt( 30.0 / ( this->getD() * t ) ) ); 

    gsl_integration_qag( &p_corr_R_F, 0.0,
			 umax,
			 err,
			 1e-8,
			 1000, GSL_INTEG_GAUSS61,
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
    gsl_sf_legendre_Pl_array( tableSize-1, cos( theta ), &PnTable[0] );

    for( Index order( 0 ); order < tableSize; ++order )
    {
	const Real Cn( 2.0 * order + 1 );

	//printf( "p %d %g\n", order, Cn * PnTable[order] * RnTable[order] );

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

    const Real p_free( this->p_free( r, r0, theta, t ) );
    const Real p_corr( this->p_corr_RnTable( theta, r, r0, t, RnTable ) ); 

    return ( p_free + p_corr ) * factor;
}

const Real PlainPairGreensFunction::drawTheta( const Real rnd,
					       const Real r, 
					       const Real r0, 
					       const Real t ) const
{
    RealVector RnTable;
    RnTable.reserve( 8 );


    {
	const unsigned int MAXORDER( 80 );
	const Real p_free_max( this->p_free( r, r0, 0.0, t ) );
	const Real integrationTolerance( p_free_max * 1e-10 );
	const Real truncationTolerance( p_free_max * 1e-8 );
    
	gsl_integration_workspace* 
	    workspace( gsl_integration_workspace_alloc( 1000 ) );
    
	Real Rn_prev( 0.0 );
	const Real RnFactor( 1.0 / ( 4.0 * M_PI * sqrt( r * r0 ) ) );
    
	unsigned int order( 0 );
	while( true ) 
	{
	    const Real Rn( this->Rn( order, r, r0, t, workspace, 
				     integrationTolerance ) );
	
	    RnTable.push_back( Rn );
	
	    //std::cerr << Rn << std::endl;
	
	    // truncate when converged enough.
	    if( fabs( Rn * RnFactor ) < truncationTolerance &&
		fabs( Rn ) < fabs( Rn_prev ) )
	    {
		break;
	    }
	
	    if( order >= MAXORDER )
	    {
		std::cerr << "Rn didn't converge." << std::endl;
		break;
	    }
	
	    Rn_prev = Rn;
	
	    ++order;
	}
    
	gsl_integration_workspace_free( workspace );
    }

    Real theta;

    {
	const Index tableSize( 200 );
	const Real thetaStep( M_PI / tableSize );
    
	//boost::array<Real, tableSize+1> pTable;
	RealVector pTable( tableSize );
    
	pTable[0] = 0.0;
	Real p_prev( 0.0 ); // p_tot_RnTable() with theta = 0 is always zero.
	Index i( 1 );
	while( true )
	{
	    Real p( this->p_tot_RnTable( thetaStep * i, r, r0, t, RnTable ) );
	    if( p < 0.0 )
	    {
		p = 0.0;
	    }
	
	    const Real value( ( p_prev + p ) * 0.5 );
	
	    pTable[i] = pTable[i-1] + value;
	
	    if( value < pTable[i] * std::numeric_limits<Real>::epsilon() ||
		i >= tableSize-1 )
	    {
		break;   // pTable is valid in [0,i].
	    }
	
	    p_prev = p;
	    ++i;
	}
    
	/*
	  const Real p_irr_r( this->p_irr( r, t, r0 ) );
	  const Real relerror( (pTable[i]*thetaStep- p_irr_r)/p_irr_r );
	  if( fabs( relerror ) >= 1e-2 )
	  {
	  std::cerr << "drawTheta: relative error estimate is large: " << relerror 
	  << "( r= " << r << ", r0= " << r0 << ", t= " << t 
	  << ")." << std::endl;
	  }
	*/
    
	const Real targetPoint( rnd * pTable[i] );
	const size_t lowerBound( gsl_interp_bsearch( &pTable[0], targetPoint, 
						     0, i ) );
	const Real low( lowerBound * thetaStep );
    
	if( pTable[lowerBound+1] - pTable[lowerBound] != 0.0 )
	{
	    theta = low + thetaStep * ( targetPoint - pTable[lowerBound] ) / 
		( pTable[lowerBound+1] - pTable[lowerBound] );
	}
	else
	{
	    // this can happen when rnd is equal to or is too close to 1.0.
	    theta = low;
	}
    }


    return theta;
}








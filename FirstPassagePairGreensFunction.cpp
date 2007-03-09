//#define NDEBUG
//#define BOOST_DISABLE_ASSERTS

#include <iostream>
#include <exception>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>

#include "bessel.hpp"

#include "HalfOrderBesselGenerator.hpp"

#include "FirstPassagePairGreensFunction.hpp"



FirstPassagePairGreensFunction::
FirstPassagePairGreensFunction( const Real D, 
				const Real kf, 
				const Real Sigma )
    :
    PairGreensFunction( D, kf, Sigma ),
    h( getkf() / ( 4.0 * M_PI * getSigma() * getSigma() * getD() ) ),
    hsigma_p_1( 1.0 + h * getSigma() ),
    a( INFINITY )
{
    this->alphaTable.reserve( 32 );
    this->expTable.reserve( 32 );
    this->psurvTable.reserve( 32 );
}

FirstPassagePairGreensFunction::~FirstPassagePairGreensFunction()
{
    ; // do nothing
}

void FirstPassagePairGreensFunction::seta( Real a )
{
    this->a = a;

    this->alphaTable.clear();
    this->expTable.clear();
}



const Real 
FirstPassagePairGreensFunction::f_alpha0( const Real alpha ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real alpha_a_m_sigma( alpha * ( a - sigma ) );
    const Real hsigma_p_1( this->hsigma_p_1 );

    Real sin_alpha_a_m_sigma;
    Real cos_alpha_a_m_sigma;
    sincos( alpha_a_m_sigma, &sin_alpha_a_m_sigma, &cos_alpha_a_m_sigma );

    const Real term1( alpha * sigma * cos_alpha_a_m_sigma );
    const Real term2( hsigma_p_1 * sin_alpha_a_m_sigma );

    const Real result( term1 + term2 );

    return result;
}

const Real 
FirstPassagePairGreensFunction::f_alpha0_aux( const Real alpha ) const

{
    const Real a( geta() );
    const Real sigma( getSigma() );

    const Real term1( ( a - sigma ) * alpha );

    const Real angle( this->hsigma_p_1 / ( sigma * alpha ) );
//    printf("%g\n",angle);
    const Real term2( std::atan( angle ) );

    const Real result( term1 - term2 );

    return result;
}

const Real 
FirstPassagePairGreensFunction::
f_alpha0_aux_df( const Real alpha ) const

{
    const Real a( geta() );
    const Real sigma( getSigma() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real alphasq( alpha * alpha );

    const Real term1( alpha - sigma );
    const Real term2( hsigma_p_1 / 
		      ( sigma * alphasq * 
			( 1.0 + ( hsigma_p_1 * hsigma_p_1 / 
				  ( sigma * sigma * alphasq ) ) ) ) );

    const Real result( term1 + term2 );

    return result;
}




const Real 
FirstPassagePairGreensFunction::
f_alpha0_aux_F( const Real alpha,
		const f_alpha0_aux_params* const params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real value( params->value );

    return gf->f_alpha0_aux( alpha ) - value;
}

const Real
FirstPassagePairGreensFunction::
f_alpha0_aux_df_F( const Real alpha,
		   const f_alpha0_aux_params* const params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 

    return gf->f_alpha0_aux_df( alpha );
}


void
FirstPassagePairGreensFunction::
f_alpha0_aux_fdf_F( const Real alpha,
		    const f_alpha0_aux_params* const params,
		    Real* const f, Real* const df )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real value( params->value );    // n * M_PI_2;

    *f = gf->f_alpha0_aux( alpha ) - value;
    *df = gf->f_alpha0_aux_df( alpha );
}




const Real 
FirstPassagePairGreensFunction::alpha0_i( const Int i ) const
{
    assert( i >= 0 );

    const Real sigma( this->getSigma() );

    const Real target( i * M_PI + M_PI_2 );
    f_alpha0_aux_params params = { this, target };


    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &f_alpha0_aux_F ),
	    &params 
	};


    // We know the range of the solution from - Pi/2 <= atan <= Pi.
    const Real rangeFactor( M_PI / ( a - sigma ) );
    Real low( i * rangeFactor );
    Real high( (i+1) * rangeFactor );

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, low, high );

    const unsigned int maxIter( 100 );

    unsigned int j( 0 );
    while( true )
    {
	gsl_root_fsolver_iterate( solver );

        low = gsl_root_fsolver_x_lower( solver );
        high = gsl_root_fsolver_x_upper( solver );
	int status( gsl_root_test_interval( low, high, 0.0, 1e-12 ) );
        //	printf("%g %g\n", low, high );


	if( status == GSL_CONTINUE )
	{
	    if( j >= maxIter )
	    {
		gsl_root_fsolver_free( solver );
		std::cerr << "alpha0_i: failed to converge." 
			  << std::endl;
		throw std::exception();
	    }
	}
	else
	{
	    break;
	}

	++j;
    }

    // printf("%d\n",j);

    const Real alpha( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );
  
    return alpha;
}



const Real 
FirstPassagePairGreensFunction::p_survival_i( const Real alpha,
					      const Real r0 ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );


    Real num1;
    {
	const Real angle_a( alpha * ( a - sigma ) );
	Real sin_a;
	Real cos_a;
	sincos( angle_a, &sin_a, &cos_a );
	num1 = alpha * sigmasq * h - 
	    alpha * ( a - sigma + a * h * sigma ) * cos_a +
	    ( hsigma_p_1 + a * sigma * alphasq ) * sin_a ;
    }

    const Real num2( num_r0( alpha, r0 ) );

    const Real den( r0 * alpha * 
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

    const Real result( num1 * num2 / den );

    return result;
}


const Real 
FirstPassagePairGreensFunction::p_leavea_i( const Real alpha,
					    const Real r0 ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    Real num1;
    {
	const Real angle_a( alpha * ( a - sigma ) );
	Real sin_a;
	Real cos_a;
	sincos( angle_a, &sin_a, &cos_a );

	num1 = - alpha * ( a - sigma + a * h * sigma ) * cos_a
	    + ( hsigma_p_1 + a * sigma * alphasq ) * sin_a;
    }
    
    Real num2;
    {
	const Real angle_r0( alpha * ( r0 - sigma ) );
	Real sin_r0;
	Real cos_r0;
	sincos( angle_r0, &sin_r0, &cos_r0 );
	
	num2 = alpha * sigma * cos_r0 + hsigma_p_1 * sin_r0;
    }
    
    
    const Real den( r0 * alpha *
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

    const Real result( num1 * num2 / den );

    return result;
}


const Real 
FirstPassagePairGreensFunction::p_leaves_i( const Real alpha,
					    const Real r0 ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );


    Real num;
    {
	const Real angle_r0( alpha * ( r0 - sigma ) );
	Real sin_r0;
	Real cos_r0;
	sincos( angle_r0, &sin_r0, &cos_r0 );
	
	num = h * sigmasq * ( alpha * sigma * cos_r0 +
			      hsigma_p_1 * sin_r0 );
    }
		      
    const Real den( r0 *
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

    const Real result( num / den );
	
    return result; // * 4 * M_PI * sigma * sigma;
}

const Real 
FirstPassagePairGreensFunction::asratio( const Real alpha,
					 const Real r0 ) const
{
    const Real a( geta() );
    const Real D( getD() );
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    const Real angle_a( alpha * ( a - sigma ) );
    Real sin_a;
    Real cos_a;
    sincos( angle_a, &sin_a, &cos_a );
    const Real num( - a * ( ( hsigma_p_1 ) * cos_a -
			    sigma * alpha * sin_a ) );
		      

    const Real den( h * sigmasq );


    const Real result( num / den );

    return result;
}


const Real
FirstPassagePairGreensFunction::num_r0( const Real alpha,
					const Real r0 ) const
{
    const Real sigma( getSigma() );
    const Real angle_r0( alpha * ( r0 - sigma ) );
    Real sin_r0;
    Real cos_r0;
    sincos( angle_r0, &sin_r0, &cos_r0 );

    const Real hsigma_p_1( this->hsigma_p_1 );
    const Real result( alpha * sigma * cos_r0 + hsigma_p_1 * sin_r0 );

    return result;
}

const Real
FirstPassagePairGreensFunction::p_int_r_i( const Real r,
					   const Real alpha,
					   const Real r0,
					   const Real num_r0 ) const
{
    // NOTE: identical to p_survival_i with a -> r.

    const Real sigma( getSigma() );

    const Real angle_r( alpha * ( r - sigma ) );
    Real sin_r;
    Real cos_r;
    sincos( angle_r, &sin_r, &cos_r );  // do sincos here; latency. 

    const Real h( geth() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    const Real num1( alpha * sigmasq * h - 
		     alpha * ( r - sigma + r * h * sigma ) * cos_r +
		     ( hsigma_p_1 + r * sigma * alphasq ) * sin_r );


//    const Real num2( num_r0( alpha, r0 ) );
    const Real num2( num_r0 );

    const Real den( r0 * alpha * 
		    ( ( r - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( r + r * h * sigma - h * sigmasq ) ) );

    const Real result( num1 * num2 / den );

    return result;
}


void 
FirstPassagePairGreensFunction::updateAlphaTable( const Real t ) const
{
    const Real a( geta() );

    RealVector& alphaTable( this->alphaTable );
    alphaTable.clear();

    alphaTable.push_back( this->alpha0_i( 0 ) );

    const Real D( getD() );
    const Real Dt2( 2.0 * D * t );
    const Real alpha0sq( alphaTable[0] * alphaTable[0] );
    const Real 
	alpha_cutoff( sqrt( gsl_sf_lambert_W0( Dt2 * alpha0sq * 
					       std::exp( Dt2 * alpha0sq ) / 
					       ( this->CUTOFF *
						 this->CUTOFF ) ) /
			    Dt2 ) );
/*    const Real alpha_cutoff( sqrt( ( - log( this->CUTOFF ) 
      / ( getD() * t ) )
      + alphaTable[0] * alphaTable[0] ) );*/


//    printf("%g %g\n", alphaTable[0], alpha_cutoff );


    const Int maxIter( 5000 );

    Int i( 1 );
    while( true )
    {
	const Real alpha_i( this->alpha0_i( i ) );
	alphaTable.push_back( alpha_i );

	if( alpha_i > alpha_cutoff )
	{
//	    printf("%d %g %g\n", 
//		   i, alpha_i, std::exp( - getD() * t * alpha_i * alpha_i ) 
//		   / alpha_i );
	    break;
	}

	if( i >= maxIter )
	{
	    std::cerr << "alphaTable: max iteration reached." << std::endl;
	    throw std::exception();
	}

	++i;
    }
}


void 
FirstPassagePairGreensFunction::updatePsurvTable( const Real r0 ) const
{
    const Real a( geta() );

    const RealVector& alphaTable( this->alphaTable );
    RealVector& psurvTable( this->psurvTable );
    psurvTable.clear();
    psurvTable.reserve( alphaTable.size() );

    for( RealVector::const_iterator i( alphaTable.begin() );
	 i != alphaTable.end(); ++i )
    {
	const Real alpha( *i );
	psurvTable.push_back( p_survival_i( alpha, r0 ) );
    }
}


void 
FirstPassagePairGreensFunction::updateNum_r0Table( RealVector& num_r0Table,
						   const Real r0 ) const
{
    const Real a( geta() );

    const RealVector& alphaTable( this->alphaTable );
    const RealVector& expTable( this->expTable );

    num_r0Table.clear();
    num_r0Table.reserve( alphaTable.size() );

    assert( alphaTable.size() >= expTable.size() );

    for( unsigned int j( 0 ); j < expTable.size(); ++j )
    {
	const Real alpha( alphaTable[j] );
	num_r0Table.push_back( num_r0( alpha, r0 ) );
    }
}


void 
FirstPassagePairGreensFunction::updateExpTable( const Real t ) const
{
    const RealVector& alphaTable( this->getAlphaTable() );

    RealVector& expTable( this->expTable );
    expTable.clear();

    const Real mDt( - getD() * t );
    const Real alpha0( alphaTable[0] );
    const Real value0( 2.0 * exp( mDt * alpha0 * alpha0 ) / alpha0 );
    expTable.push_back( value0 );

    const Real factor( 1.0 / value0 );

    for( RealVector::size_type j( 1 ); j < alphaTable.size(); ++j )
    {
	const Real alpha( alphaTable[j] );
	const Real value( 2.0 * std::exp( mDt * alpha * alpha ) / alpha );
	expTable.push_back( value );

	if( value * factor < this->CUTOFF )
	{
	    // printf("cutoff at %d; %g\n", j, value );
	    break;
	}
    }

}


const Real 
FirstPassagePairGreensFunction::p_survival( const Real t,
					    const Real r0 ) const
{
    Real p( 0.0 );

    const RealVector& alphaTable( this->getAlphaTable() );
    const RealVector& expTable( this->expTable );
    const RealVector& psurvTable( this->psurvTable );

    assert( alphaTable.size() >= expTable.size() );

    for( RealVector::size_type j( 0 ); j < expTable.size(); ++j )
    {
	const Real value( psurvTable[j] );
	p += value * expTable[j];
    }

    return p;
}

const Real 
FirstPassagePairGreensFunction::p_leaves( const Real t,
					  const Real r0 ) const
{
    Real p( 0.0 );

    const RealVector& alphaTable( this->getAlphaTable() );
    const RealVector& expTable( this->expTable );

    this->updateExpTable( t );

    for( unsigned int j( 0 ); j < expTable.size(); ++j )
    {
	const Real value( p_leaves_i( alphaTable[j], r0 ) );

	p += value * expTable[j];
    }

    return p;
}

const Real 
FirstPassagePairGreensFunction::p_leavea( const Real t,
					  const Real r0 ) const
{
    Real p( 0.0 );

    const RealVector& alphaTable( this->getAlphaTable() );
    const RealVector& expTable( this->expTable );

    this->updateExpTable( t );

    for( unsigned int j( 0 ); j < expTable.size(); ++j )
    {
	const Real value( p_leavea_i( alphaTable[j], r0 ) );

	p += value * expTable[j];
    }

    return p;
}



const Real 
FirstPassagePairGreensFunction::p_int_r( const Real r,
					 const Real t,
					 const Real r0,
					 const RealVector& num_r0Table ) const
{
    Real p( 0.0 );

    const RealVector& alphaTable( this->getAlphaTable() );
    const RealVector& expTable( this->expTable );

    assert( alphaTable.size() >= expTable.size() );

    for( RealVector::size_type j( 0 ); j < expTable.size(); ++j )
    {
	const Real alpha( alphaTable[j] );
	const Real value( p_int_r_i( r, alpha, r0, num_r0Table[j] ) );
	p += value * expTable[j];
    }

    return p;
}





const Real
FirstPassagePairGreensFunction::p_survival_F( const Real t,
					      const p_survival_params* params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real r0( params->r0 );
    const Real rnd( params->rnd );

    gf->updateExpTable( t );

    return gf->p_survival( t, r0 ) - rnd;
}


const Real
FirstPassagePairGreensFunction::p_int_r_F( const Real r,
					   const p_int_r_params* params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real t( params->t );
    const Real r0( params->r0 );
    const Real psurv( params->psurv );
    const RealVector& num_r0Table( params->num_r0Table );
    const Real rnd( params->rnd );

    return ( gf->p_int_r( r, t, r0, num_r0Table ) / psurv ) - rnd;
}


const Real FirstPassagePairGreensFunction::drawTime( const Real rnd, 
						     const Real r0 ) const
{
    assert( rnd <= 1.0 && rnd >= 0.0 );
    assert( r0 > 0.0 );

    p_survival_params params = { this, r0, rnd };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &p_survival_F ),
	    &params 
	};

    Real low( 1e-5 );
    Real high( 10.0 );

    this->updateAlphaTable( low );
    this->updatePsurvTable( r0 );

    //FIXME: adjust high here.

    // adjust low to make sure tha f( low ) and f( high ) straddle.
    const Real highvalue( GSL_FN_EVAL( &F, high ) );
    while( GSL_FN_EVAL( &F, low ) * highvalue >= 0.0 )
    {
	low *= .1;
	//printf("drawTime: adjusting low: %g\n",low);

	if( fabs( low ) <= 1e-50 )
	{
	    std::cerr << "Couldn't adjust low. (" << low <<
		")" << std::endl;
	    throw std::exception();
	    
	}
	this->updateAlphaTable( low );
	this->updatePsurvTable( r0 );
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
	int status( gsl_root_test_interval( low, high, .0, this->CUTOFF ) );

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

#ifndef NDEBUG
    Real p( p_survival( t, r0 ) );
    Real leaves( p_leaves( t, r0 ) );
    Real leavea( p_leavea( t, r0 ) );
    const Real error( p - leaves - leavea );
    if( fabs( error ) >= this->CUTOFF )
    {
	printf("ERROR >= CUTOFF: %g; %g %g %g\n", error, p, leaves, leavea );
    }
#endif

    return t;
}

const bool FirstPassagePairGreensFunction::drawEventType( const Real rnd, 
							  const Real r0,
							  const Real t ) const
{
    // NOTE: The huge assumption for this method to operate correctly is that
    // drawTime() was called immediately before invokation of this with
    // the same parameter r0.

    Real p( p_survival( t, r0 ) );
    Real leaves( p_leaves( t, r0 ) );

    Real value( leaves / p );

//    Real leavea( p - leaves );
//    printf("%g %g %g %g\n", value, p, leaves, leavea);

    if( rnd < value )  
    {
	return true;   // leaves
    }
    else 
    {
	return false;  // leavea
    }
}


const Real FirstPassagePairGreensFunction::drawR( const Real rnd, 
						  const Real r0, 
						  const Real t ) const
{
    assert( rnd <= 1.0 && rnd >= 0.0 );
    assert( r0 > 0.0 );

    this->updateAlphaTable( t );
    this->updateExpTable( t );
    this->updatePsurvTable( r0 );

    const Real psurv( p_survival( t, r0 ) );

    RealVector num_r0Table;
    updateNum_r0Table( num_r0Table, r0 );

    p_int_r_params params = { this, t, r0, psurv, num_r0Table, rnd };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &p_int_r_F ),
	    &params 
	};

    Real low( this->getSigma() );
    Real high( this->geta() );

    //const Real lowvalue( GSL_FN_EVAL( &F, low*1.1 ) );
    //    const Real highvalue( GSL_FN_EVAL( &F, high*0.9 ) );
    //printf("%g %g %g\n", lowvalue, highvalue, psurv );

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
	int status( gsl_root_test_interval( low, high, .0, this->CUTOFF ) );

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


    Real r( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );

    return r;
}


const Real FirstPassagePairGreensFunction::f_alpha( const Real alpha,
						    const Int n ) const
{
    const Real a( this->geta() );
    const Real aAlpha( a * alpha );
    const Real sigmaAlpha( getSigma() * alpha );
    const Real hSigma( geth() * getSigma() );
    const Real realn( static_cast<Real>( n ) );

    const Real hSigma_m_n( hSigma - realn );

#if 0
    // Numerical recipes
    Real tmp, jas1, yas1, jas2, yas2, jaa, yaa;
    bessjy( sigmaAlpha, realn + 0.5, &jas1, &yas1, &tmp, &tmp );
    bessjy( sigmaAlpha, realn + 1.5, &jas2, &yas2, &tmp, &tmp );
    bessjy( aAlpha, realn + 0.5, &jaa, &yaa, &tmp, &tmp );
#else
    // GSL
    const Real factors( sqrt( sigmaAlpha * M_2_PI ) );
    const Real factora( sqrt( aAlpha * M_2_PI ) );
    const Real jas1( gsl_sf_bessel_jl( n, sigmaAlpha ) * factors );
    const Real yas1( gsl_sf_bessel_yl( n, sigmaAlpha ) * factors );
    const Real jas2( gsl_sf_bessel_jl( n+1, sigmaAlpha ) * factors );
    const Real yas2( gsl_sf_bessel_yl( n+1, sigmaAlpha ) * factors );
    const Real jaa( gsl_sf_bessel_jl( n, aAlpha ) * factora );
    const Real yaa( gsl_sf_bessel_yl( n, aAlpha ) * factora );
#endif

//    printf("b %g %g %g %g %g %g\n", jas1, jas2, yas1, yas2, yaa, jaa );
    const Real term1( ( hSigma_m_n * jas1 + sigmaAlpha * jas2 ) * yaa );
    const Real term2( ( hSigma_m_n * yas1 + sigmaAlpha * yas2 ) * jaa );

//    printf("s %g %g %g %g\n", hSigma_m_n * jas1 * yaa, sigmaAlpha * jas2 * yaa,
//	   hSigma_m_n * yas1 * jaa, sigmaAlpha * yas2 * jaa);

//    printf("t %g %g %g %g\n", alpha, term1, term2, term1-term2 );// cos(f_alpha_aux( alpha,n )) );
    const Real result( term1 - term2 );
    
    return result;
}

const Real G( const unsigned int n, const unsigned int k )
{
    //std::cerr << n << ' ' << k << std::endl;
    return gsl_sf_fact( n + k ) / ( gsl_sf_fact( k ) * gsl_sf_fact( n - k ) );
}

const Real FirstPassagePairGreensFunction::P( const Int n,
					      const Real x )
{
    Real result( 0.0 );

    Real sx2( 1.0 );
    const Real x2sq_r( 1.0 / gsl_pow_2( x + x ) );

    const unsigned int maxm( n/2 );
    for( unsigned int m( 0 ); m <= maxm; ++m )
    {
	const Real term1( 1.0 + (-2.0) * ( m % 2 ) ); // (-1)^m
	const Real value( term1 * sx2 * G( n, 2 * m ) );
	result += value;

	sx2 *= x2sq_r;
    }

    return result;
}


const Real FirstPassagePairGreensFunction::Q( const Int n,
					      const Real x )
{
    Real result( 0.0 );

    Real sx2( 1.0 / ( x + x ) );
    const Real x2sq( sx2 * sx2 );

    const unsigned int maxm( (n+1)/2 ); // sum_(0)^((n-1)/2)
    for( unsigned int m( 0 ); m < maxm; ++m )
    {
	const Real term1( 1.0 + (-2.0) * ( m % 2 ) ); // (-1)^m

	const Real value( term1 * sx2 * G( n, 2 * m + 1 ) );
	result += value;

	sx2 *= x2sq;
    }

    return result;
}

const Real 
FirstPassagePairGreensFunction::f_alpha_aux( const Real alpha, 
					     const Int n ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );

    const Real aAlpha( a * alpha );
    const Real sigmaAlpha( sigma * alpha );

    const Real realn( static_cast<Real>( n ) );
    const Real n_m_hSigma( n - h * sigma );

    /*(a - s) u - 
      ArcTan[( P[n, a u] ((-n + h s) P[n, s u] + s u Q[1 + n, s u]) -
               Q[n, a u] (s u P[1 + n, s u] + (n - h s) Q[n, s u]) )/
             ( Q[n, a u] ((-n + h s) P[n, s u] + s u Q[1 + n, s u]) + 
               P[n, a u] (s u P[1 + n, s u] + (n - h s) Q[n, s u]) )]
    */

    const Real term1( ( a - sigma ) * alpha );

    const Real Pa( P( n, aAlpha ) );
    const Real Ps( P( n, sigmaAlpha ) );
    const Real Psp( P( n+1, sigmaAlpha ) );
    const Real Qa( Q( n, aAlpha ) );
    const Real Qs( Q( n, sigmaAlpha ) );
    const Real Qsp( Q( n+1, sigmaAlpha ) );

    //printf("%g %g %g %g %g %g\n",Pa,Ps,Qa,Qs,Psp,Qsp);

    /*
    const Real sigmaAlphaQs_m_n_m_hSigmaPs( sigmaAlpha * Qs - 
					      n_m_hSigma * Ps );
    const Real n_m_hSigmaQs_p_sigmaAlphaPs( n_m_hSigma * Qs + 
					      sigmaAlpha * Ps );
    */

    const Real PaQa( Pa / Qa );

    const Real t1( Qa * ( sigmaAlpha * Qsp - n_m_hSigma * Ps ) ); 
    const Real t2( Pa * ( n_m_hSigma * Qs  + sigmaAlpha * Psp ) );

    const Real t3( Qa * sigmaAlpha * Qsp + Pa * sigmaAlpha * Psp );
    const Real t4( Pa * n_m_hSigma * Qs  - Qa * n_m_hSigma * Ps );
//    printf("ppp %f %f %g %g %g %g\n", Qa * sigmaAlpha * Qsp,Pa * sigmaAlpha * Psp,Pa * n_m_hSigma * Qs, Qa * n_m_hSigma * Ps );
//    printf("ppp %f %f %g %g %g %g\n", t1, t2, t1 + t2, t3, t4, t3+t4 );


    const Real angle( ( ( PaQa * sigmaAlpha * Qsp - PaQa * n_m_hSigma * Ps ) - 
			( n_m_hSigma * Qs  + sigmaAlpha * Psp ) ) /
		      ( ( sigmaAlpha * Qsp - n_m_hSigma * Ps ) + 
			( PaQa * n_m_hSigma * Qs  + PaQa * sigmaAlpha * Psp ) ) );
//    printf("%g\n", angle );

/*
    const Real angle( ( Pa * ( sigmaAlpha * Qsp - n_m_hSigma * Ps ) - 
			Qa * ( n_m_hSigma * Qs  + sigmaAlpha * Psp ) ) /
		      ( Qa * ( sigmaAlpha * Qsp - n_m_hSigma * Ps ) + 
			Pa * ( n_m_hSigma * Qs  + sigmaAlpha * Psp ) ) );
*/
    const Real term2( std::atan( angle ) );


/*
    const Real term2( std::atan2( 
			  ( Pa * ( sigmaAlpha * Qsp - n_m_hSigma * Ps ) - 
			    Qa * ( n_m_hSigma * Qs  + sigmaAlpha * Psp ) ),
			  ( Qa * ( sigmaAlpha * Qsp - n_m_hSigma * Ps ) + 
			    Pa * ( n_m_hSigma * Qs  + sigmaAlpha * Psp ) ) ) );
*/
    const Real result( term1 - term2 );
//    printf("value %g cos %g result %g \n",f_alpha(alpha,n), cos( result ), result);
//    printf("aux %18.18g %g %g %g %g\n", alpha, result, term1, term2, angle );
    

    return result;
}


const Real 
FirstPassagePairGreensFunction::
f_alpha_aux_F( const Real alpha,
	       const f_alpha_aux_params* const params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Int n( params->n );
    const Real value( params->value );

    return gf->f_alpha_aux( alpha, n ) - value;
    //return gf->f_alpha( alpha, n );// - value;
}


const Real 
FirstPassagePairGreensFunction::alpha_i( const Int i, const Int n ) const
{
    assert( i >= 0 );

    const Real sigma( this->getSigma() );

    Real target( (i+1) * M_PI + M_PI_2 ); //+ n/2 * M_PI );
    f_alpha_aux_params params = { this, n, target };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &f_alpha_aux_F ),
	    &params 
	};


    // We know the range of the solution from - Pi/2 <= atan <= Pi.
    const Real alphaMid( target / ( a - sigma ) );
    const Real alphaHalfRange( M_PI_2 / ( a - sigma ) );
    Real low( alphaMid - alphaHalfRange );
    Real high( alphaMid + alphaHalfRange );
    printf("target %g low %g high %g\n",target, low, high);

    Real lowvalue;
    Real highvalue;

    while( true )
    {
	Real lowvalue( f_alpha(low,n) );
	Real highvalue( f_alpha(high,n) );

	if( lowvalue * highvalue >= 0 )
	{
	    printf("lh: %g %g %g %g\n", low, high, lowvalue, highvalue );
	    target += M_PI;
	    low += (target - M_PI_2) / (a-sigma);
	    high += (target + M_PI_2) / (a-sigma);
	}
	else
	{
	    printf("ok: %g %g %g %g\n", low, high, lowvalue, highvalue );
	    break;
	}
    }



    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, low, high );

    const unsigned int maxIter( 100 );

    unsigned int j( 0 );
    while( true )
    {
	gsl_root_fsolver_iterate( solver );

        low = gsl_root_fsolver_x_lower( solver );
        high = gsl_root_fsolver_x_upper( solver );
	int status( gsl_root_test_interval( low, high, 0.0, 1e-15 ) );

	if( status == GSL_CONTINUE )
	{
	    if( j >= maxIter )
	    {
		gsl_root_fsolver_free( solver );
		std::cerr << "alpha_i: failed to converge." 
			  << std::endl;
		throw std::exception();
	    }
	}
	else
	{
	    break;
	}

	++j;
    }


    const Real alpha( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );
  
    printf("%d %g\n",j,alpha);

    return alpha;
}


    
const Real FirstPassagePairGreensFunction::drawTheta( const Real rnd,
						      const Real r, 
						      const Real r0, 
						      const Real t ) const
{

    for( int i(0); i< 10; ++i )
    {
	printf("%d %g %g\n", i, P(i,.2), Q(i,.2) );

    }


    return 0.0;
}
    

//#define NDEBUG
//#define BOOST_DISABLE_ASSERTS

#include <iostream>
#include <exception>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_lambert.h>
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
FirstPassagePairGreensFunction::f_alpha_survival( const Real alpha ) const
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
FirstPassagePairGreensFunction::f_alpha_survival_aux( const Real alpha ) const

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
f_alpha_survival_aux_df( const Real alpha ) const

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
f_alpha_survival_aux_F( const Real alpha,
			const f_alpha_survival_aux_params* const params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real value( params->value );

    return gf->f_alpha_survival_aux( alpha ) - value;
}

const Real
FirstPassagePairGreensFunction::
f_alpha_survival_aux_df_F( const Real alpha,
			   const f_alpha_survival_aux_params* const params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 

    return gf->f_alpha_survival_aux_df( alpha );
}


void
FirstPassagePairGreensFunction::
f_alpha_survival_aux_fdf_F( const Real alpha,
			    const f_alpha_survival_aux_params* const params,
			    Real* const f, Real* const df )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real value( params->value );    // n * M_PI_2;

    *f = gf->f_alpha_survival_aux( alpha ) - value;
    *df = gf->f_alpha_survival_aux_df( alpha );
}




const Real 
FirstPassagePairGreensFunction::alpha_survival_n( const Int n ) const
{
    assert( n >= 0 );

    const Real sigma( this->getSigma() );

    const Real target( n * M_PI + M_PI_2 );
    f_alpha_survival_aux_params params = { this, target };


    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &f_alpha_survival_aux_F ),
	    &params 
	};


    // We know the range of the solution from - Pi/2 <= atan <= Pi.
    const Real rangeFactor( M_PI / ( a - sigma ) );
    Real low( n * rangeFactor );
    Real high( low + rangeFactor );

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, low, high );

    const unsigned int maxIter( 1000 );

    unsigned int i( 0 );
    while( true )
    {
	gsl_root_fsolver_iterate( solver );

        low = gsl_root_fsolver_x_lower( solver );
        high = gsl_root_fsolver_x_upper( solver );
	int status( gsl_root_test_interval( low, high, 0.0, 1e-12 ) );
        //	printf("%g %g\n", low, high );


	if( status == GSL_CONTINUE )
	{
	    if( i >= maxIter )
	    {
		gsl_root_fsolver_free( solver );
		std::cerr << "alpha_survival_n: failed to converge." 
			  << std::endl;
		throw std::exception();
	    }
	}
	else
	{
	    break;
	}

	++i;
    }

    // printf("%d\n",i);

    Real alpha( gsl_root_fsolver_root( solver ) );
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


    Real num2;
    {
	const Real angle_r0( alpha * ( r0 - sigma ) );
	Real sin_r0;
	Real cos_r0;
	sincos( angle_r0, &sin_r0, &cos_r0 );
	num2 = alpha * sigma * cos_r0 + hsigma_p_1 * sin_r0 ;
    }
		      

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




void 
FirstPassagePairGreensFunction::updateAlphaTable( const Real t ) const
{
    const Real a( geta() );

    RealVector& alphaTable( this->alphaTable );
    alphaTable.clear();

    alphaTable.push_back( this->alpha_survival_n( 0 ) );

    const Real D( getD() );
    const Real Dt2( 2.0 * D * t );
    const Real alpha0sq( alphaTable[0] * alphaTable[0] );
    const Real 
	alpha_cutoff( sqrt( gsl_sf_lambert_W0( Dt2 * alpha0sq * 
					       std::exp( Dt2 * alpha0sq ) / 
					       ( this->CUTOFF *
						 this->CUTOFF ) ) /
			    Dt2 ) );
/*    const Real alpha_cutoff( sqrt( ( - log( this->ALPHA_CUTOFF ) 
				     / ( getD() * t ) )
				     + alphaTable[0] * alphaTable[0] ) );*/


//    printf("%g %g\n", alphaTable[0], alpha_cutoff );


    const Int maxIter( 5000 );

    Int i( 1 );
    while( true )
    {
	const Real alpha_i( this->alpha_survival_n( i ) );
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
FirstPassagePairGreensFunction::updateExpTable( const Real t ) const
{
    const RealVector& alphaTable( this->getAlphaTable() );

    RealVector& expTable( this->expTable );
    expTable.clear();
    expTable.reserve( 8 );

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

    this->updateExpTable( t );

    for( unsigned int j( 0 ); j < expTable.size(); ++j )
    {
	const Real value( p_survival_i( alphaTable[j], r0 ) );

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




const Real FirstPassagePairGreensFunction::f_alpha( const Real alpha,
						    const Int n ) const
{
    const Real a( this->geta() );
    const Real aalpha( a * alpha );
    const Real alphaSigma( alpha * getSigma() );
    const Real hSigma( geth() * getSigma() );
    const Real realn( static_cast<Real>( n ) );

    const Real hSigma_m_n( hSigma - realn );

    Real tmp, jas1, yas1, jas2, yas2, jaa, yaa;

    bessjy( alphaSigma, realn + 0.5, &jas1, &yas1, &tmp, &tmp );
    bessjy( alphaSigma, realn + 1.5, &jas2, &yas2, &tmp, &tmp );
    bessjy( aalpha, realn + 0.5, &jaa, &yaa, &tmp, &tmp );


    const Real term1( ( hSigma_m_n * jas1 + alphaSigma * jas2 ) * yaa );
    const Real term2( ( hSigma_m_n * yas1 + alphaSigma * yas2 ) * jaa );

    const Real result( term1 - term2 );
    
    return result;
}


const Real
FirstPassagePairGreensFunction::p_survival_F( const Real t,
					      const p_survival_params* params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real r0( params->r0 );
    const Real rnd( params->rnd );

    return gf->p_survival( t, r0 ) - rnd;
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
	this->updateExpTable( low );
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
  
    Real t( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );

    Real p = p_survival( t, r0 );
    Real leaves( p_leaves( t, r0 ) );
    Real leavea( p_leavea( t, r0 ) );

#ifndef NDEBUG
    const Real error( p - leaves - leavea );
    if( fabs( error ) >= this->CUTOFF )
    {
	printf("ERROR >= CUTOFF: %g; %g %g %g\n", error, p, leaves, leavea );
    }
#endif

    return t;
}

const Real FirstPassagePairGreensFunction::drawR( const Real rnd, 
						  const Real r0, 
						  const Real t ) const
{


    return 0.0;
}
    
const Real FirstPassagePairGreensFunction::drawTheta( const Real rnd,
						      const Real r, 
						      const Real r0, 
						      const Real t ) const
{


    return 0.0;
}
    

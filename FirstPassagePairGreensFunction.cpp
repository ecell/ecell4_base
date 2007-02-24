#define HAVE_INLINE

//#define NDEBUG
//#define BOOST_DISABLE_ASSERTS

#include <exception>
#include <vector>

#include <boost/array.hpp>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
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
    hsigma_p_1( 1.0 + h * getSigma() )
{
    ; // do nothing
}

FirstPassagePairGreensFunction::~FirstPassagePairGreensFunction()
{
    ; // do nothing
}


const Real 
FirstPassagePairGreensFunction::f_alpha_survival( const Real alpha,
						  const Real a ) const
{
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real alpha_a_m_sigma( alpha * ( a - sigma ) );

    const Real term1( alpha * sigma * cos( alpha_a_m_sigma ) );
    const Real term2( ( 1.0 + h * sigma ) * sin( alpha_a_m_sigma ) );

    const Real result( term1 + term2 );

    return result;
}

const Real 
FirstPassagePairGreensFunction::f_alpha_survival_aux( const Real alpha,
						      const Real a ) const

{
    const Real sigma( getSigma() );

    const Real term1( ( a - sigma ) * alpha );
    const Real term2( atan( this->hsigma_p_1 / ( sigma * alpha ) ) );

    const Real result( term1 - term2 );

    return result;
}

const Real 
FirstPassagePairGreensFunction::f_alpha_survival_aux_df( const Real alpha,
							 const Real a ) const

{
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
    const Real a( params->a );
    const Real value( params->value );    // n * M_PI_2;

    return gf->f_alpha_survival_aux( alpha, a ) - value;
}

const Real
FirstPassagePairGreensFunction::
f_alpha_survival_aux_df_F( const Real alpha,
			   const f_alpha_survival_aux_params* const params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real a( params->a );
    //    const Real value( params->value );    // n * M_PI_2;

    return gf->f_alpha_survival_aux_df( alpha, a );
}


void
FirstPassagePairGreensFunction::
f_alpha_survival_aux_fdf_F( const Real alpha,
			    const f_alpha_survival_aux_params* const params,
			    Real* const f, Real* const df )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real a( params->a );
    const Real value( params->value );    // n * M_PI_2;

    *f = gf->f_alpha_survival_aux( alpha, a ) - value;
    *df = gf->f_alpha_survival_aux_df( alpha, a );
}




const Real 
FirstPassagePairGreensFunction::alpha_survival_n( const Real a,
						  const Int n ) const
{
    assert( a > 0 && n >= 0 );

    const Real sigma( this->getSigma() );

    const Real target( n * M_PI + M_PI_2 );
    f_alpha_survival_aux_params params = { this, a, target };


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
//    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_bisection );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, low, high );

    const unsigned int maxIter( 1000 );

    unsigned int i( 0 );
    Real alpha;
    while( true )
    {
	gsl_root_fsolver_iterate( solver );

	alpha = gsl_root_fsolver_root( solver );
        low = gsl_root_fsolver_x_lower( solver );
        high = gsl_root_fsolver_x_upper( solver );
	int status( gsl_root_test_interval( low, high, 0.0, 1e-15 ) );
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

    gsl_root_fsolver_free( solver );
  
    return alpha;
}



const Real 
FirstPassagePairGreensFunction::p_survival_i( const Real alpha,
					      const Real t,
					      const Real r0,
					      const Real a ) const
{
    const Real sigma( getSigma() );
    const Real h( geth() );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    const Real hs_p_1( h * sigma + 1.0 );

//    const Real term1( 2.0 * exp( - D * t * alphasq ) );

    const Real angle_a( alpha * ( a - sigma ) );

    const Real term2( alpha * sigmasq * h - 
		      alpha * ( a - sigma + a * h * sigma ) * cos( angle_a ) +
		      ( hs_p_1 + a * sigma * alphasq ) * sin( angle_a ) );


    const Real angle_r0( alpha * ( r0 - sigma ) );
    const Real term3( alpha * sigma * cos( angle_r0 ) +
		      hs_p_1 * sin( angle_r0 ) );
		      

    const Real den( r0 * alphasq * 
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hs_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

    const Real result( term2 * term3 / den );

    return result;
}


const Real 
FirstPassagePairGreensFunction::p_leavea_i( const Real alpha,
					    const Real t,
					    const Real r0,
					    const Real a ) const
{
    const Real sigma( getSigma() );
    const Real h( geth() );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    const Real hs_p_1( h * sigma + 1.0 );

//    const Real term1( 2.0 * exp( - D * t * alphasq ) );

    const Real angle_a( alpha * ( a - sigma ) );

    const Real term2( - alpha * ( a - sigma + a * h * sigma ) * cos( angle_a ) 
		      + ( hs_p_1 + a * sigma * alphasq ) * sin( angle_a ) );

    const Real angle_r0( alpha * ( r0 - sigma ) );
    const Real term3( alpha * sigma * cos( angle_r0 ) +
		      hs_p_1 * sin( angle_r0 ) );
		      

    const Real den( r0 * alphasq *
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hs_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

    const Real result( term2 * term3 / den );


    return result;// * 4 * M_PI * a * a;
}


const Real 
FirstPassagePairGreensFunction::p_leaves_i( const Real alpha,
					    const Real t,
					    const Real r0,
					    const Real a ) const
{
    const Real sigma( getSigma() );
    const Real h( geth() );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    const Real hs_p_1( h * sigma + 1.0 );


//    const Real term1( 2.0 * exp( - D * t * alphasq ) );

    const Real angle_r0( alpha * ( r0 - sigma ) );
    const Real term3( h * sigmasq * ( alpha * sigma * cos( angle_r0 ) +
				      hs_p_1 * sin( angle_r0 ) ) );
		      

    const Real den( r0 * alpha *
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hs_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

    const Real result( term3 / den );
	
    return result; // * 4 * M_PI * sigma * sigma;
}

const Real 
FirstPassagePairGreensFunction::asratio( const Real alpha,
					 const Real t,
					 const Real r0,
					 const Real a ) const
{
    const Real D( getD() );
    const Real sigma( getSigma() );
    const Real h( geth() );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    const Real hs_p_1( h * sigma + 1.0 );

//    const Real term1( 2.0 * exp( - D * t * alphasq ) );

    const Real angle_a( alpha * ( a - sigma ) );
    const Real term3( - a * ( ( hs_p_1 ) * cos( angle_a ) -
			    sigma * alpha * sin( angle_a ) ) );
		      

    const Real den( h * sigmasq );


    const Real result( term3 / den );

    return result;
}




void 
FirstPassagePairGreensFunction::updateAlphaTable( RealVector& alphaTable,
						  const Real t,
						  const Real a ) const
{
    if( alphaTable.empty() )
    {
	alphaTable.push_back( this->alpha_survival_n( a, 0 ) );
    }

    const Real D( getD() );
    const Real Dt2( 2.0 * D * t );
    const Real alpha0sq( alphaTable[0] * alphaTable[0] );
    const Real 
	alpha_cutoff( sqrt( gsl_sf_lambert_W0( Dt2 * alpha0sq * 
					       exp( Dt2 * alpha0sq ) / 
					       ( this->ALPHA_CUTOFF *
						 this->ALPHA_CUTOFF ) ) /
			    Dt2 ) );
/*    const Real alpha_cutoff( sqrt( ( - log( this->ALPHA_CUTOFF ) 
				     / ( getD() * t ) )
				     + alphaTable[0] * alphaTable[0] ) );*/


    printf("%g %g\n", alphaTable[0], alpha_cutoff );


    const Int maxIter( 1000 );

    Int i( 1 );
    while( true )
    {
	const Real alpha_i( this->alpha_survival_n( a, i ) );
	alphaTable.push_back( alpha_i );

	printf("%d %g %g\n", 
	       i, alpha_i, exp( - getD() * t * alpha_i * alpha_i ) / alpha_i );

	if( alpha_i > alpha_cutoff )
	{
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

const Real 
FirstPassagePairGreensFunction::p_survival( const Real t,
					    const Real r0,
					    const Real a
					    ) const
{
    RealVector alphaTable;
    alphaTable.reserve( 8 );

    this->updateAlphaTable( alphaTable, t, a );
    
    const Real D( this->getD() );

    RealVector expTable( alphaTable.size() );
    for( RealVector::size_type j( 0 ); j != expTable.size(); ++j )
    {
	const Real alphasq( alphaTable[j] * alphaTable[j] );
	expTable[j] = 2.0 * exp( - D * t * alphasq );
    }


    Real p( 0.0 );
    Real p_a( 0.0 );
    Real p_s( 0.0 );
    for( unsigned int j( 0 ); j < alphaTable.size(); ++j )
    {
	const Real value( p_survival_i( alphaTable[j], t, r0, a ) );
	const Real value_a( p_leavea_i( alphaTable[j], t, r0, a ) );
	const Real value_s( p_leaves_i( alphaTable[j], t, r0, a ) );

	p += value * expTable[j];
	p_a += value_a * expTable[j];
	p_s += value_s * expTable[j];

	printf("%d %g %g %g %g %g\n", j, value, p, p_a + p_s, p_a, p_s );
    }

    return p;
}




const Real FirstPassagePairGreensFunction::f_alpha( const Real alpha,
						    const Real a,
						    const Int n ) const
{
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



const Real FirstPassagePairGreensFunction::drawTime( const Real rnd, 
						     const Real r0,
						     const Real maxt ) const
{

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
    

#define HAVE_INLINE

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

#include "HalfOrderBesselGenerator.hpp"

#include "FirstPassagePairGreensFunction.hpp"



FirstPassagePairGreensFunction::
FirstPassagePairGreensFunction( const Real D, 
				const Real kf, 
				const Real Sigma )
    :
    PairGreensFunction( D, kf, Sigma ),
    h( getkf() / ( 4.0 * M_PI * getSigma() * getSigma() * getD() ) )
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
    const Real h( geth() );
    const Real hsigma_p_1( 1.0 + h * sigma );

    const Real term1( ( a - sigma ) * alpha );
    const Real term2( atan( hsigma_p_1 / ( sigma * alpha ) ) );
//    const Real term2( atan2( hsigma_p_1, ( sigma * alpha ) ) );
    
    const Real result( term1 - term2 );

    return result;
}

const Real 
FirstPassagePairGreensFunction::f_alpha_survival_aux_df( const Real alpha,
							 const Real a ) const

{
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real hsigma_p_1( 1.0 + h * sigma );

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
						  const Int n,
						  const Real lower ) const
{
    assert( lower > 0 );

    const Real target( n * M_PI + M_PI_2 );
    f_alpha_survival_aux_params params = { this, a, target };


    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &f_alpha_survival_aux_F ),
	    &params 
	};

    Real factor( 2.5 );

    Real low( lower );
    Real high( lower * factor );

    // adjust low to make sure tha f( low ) and f( high ) straddle.
    const Real lowvalue( GSL_FN_EVAL( &F, low ) );
    while( GSL_FN_EVAL( &F, high ) * lowvalue >= 0.0 )
    {
	printf("alpha_survival_n: adjusting high: %g\n",high);
	factor *= factor;
	high *= factor;
	if( fabs( low ) <= 1e-50 )
	{
	    std::cerr << "Couldn't adjust high. (" << low <<
		      ")" << std::endl;
	    throw std::exception();
	    
	}
    }


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

    printf("%d\n",i);

    gsl_root_fsolver_free( solver );
  
    return alpha;
}



const Real 
FirstPassagePairGreensFunction::p_survival_i( const Real r,
					      const Real t,
					      const Real r0,
					      const Real alpha,
					      const Real a ) const
{
    const Real D( getD() );
    const Real sigma( getSigma() );
    const Real h( geth() );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    const Real alpha_a_m_sigma( alpha * ( a - sigma ) );

    const Real hs_p_1( h * sigma + 1.0 );

    const Real term1( exp( - D * t * alphasq ) );
    const Real term2( alpha * sigma * cos( alpha * ( r - sigma ) ) +
		      hs_p_1 * cos( alpha * ( r - sigma ) ) );
    const Real term3( alpha * sigma * cos( alpha * ( r0 - sigma ) ) +
		      hs_p_1 * cos( alpha * ( r0 - sigma ) ) );
		      

    const Real PI2rr0( 2.0 * M_PI * r * r0 );

    const Real den( PI2rr0 * 
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hs_p_1 * ( - h * sigmasq + a * h * sigma + a ) ) );


    const Real result( term1 * term2 * term3 / den );

    return result;
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
    

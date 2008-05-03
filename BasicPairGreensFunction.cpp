//#define NDEBUG
//#define BOOST_DISABLE_ASSERTS

#include <exception>
#include <vector>
#include <sstream>
#include <cmath>

#include <boost/bind.hpp>


#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
//#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_roots.h>

//#include "bessel.hpp"

#include "freeFunctions.hpp"

#include "funcSum.hpp"

//#include "HalfOrderBesselGenerator.hpp"

#include "BasicPairGreensFunction.hpp"



BasicPairGreensFunction::BasicPairGreensFunction( const Real D, 
						  const Real kf, 
						  const Real Sigma )
    :
    PairGreensFunction( D, kf, Sigma ),
    kD( 4.0 * M_PI * getSigma() * getD() ),
    alpha( ( 1.0 + ( getkf() / getkD() ) ) * ( sqrt( getD() ) / getSigma() ) )
{
    ; // do nothing
}

BasicPairGreensFunction::~BasicPairGreensFunction()
{
    ; // do nothing
}




#if 0 
const Real 
BasicPairGreensFunction::p_corr_R2( const Real alpha, 
				    const p_corr_R2_params* const params )
{
    const Real SIGMA2KFp1( 2.0 * params->Sigma * params->kf + 1.0 );
    const Real SIGMA2U( 2.0 * alpha * params->Sigma );

    const Real costheta( cos( params->theta ) );
    
    const int order_min( params->order );
    const int order_max( params->order_max );
    HalfOrderBesselGenerator us( alpha * params->Sigma, order_min, order_max );
    HalfOrderBesselGenerator ur( alpha * params->r, order_min, order_max );
    HalfOrderBesselGenerator ur0( alpha * params->r0, order_min, order_max );
    
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

    const Real exp_term( exp( - params->D * alpha * alpha * params->t ) * alpha );
    result *= exp_term;



    
    if( isnan( result ) )
    {
	std::cerr << "NaN in p_corr_R" << std::endl;
	std::cerr << alpha << ' ' << order_min << ' ' << exp_term << std::endl;
	//      std::cerr << "R1F1 " << R1F1 << " R2 F2" << R2 << ' ' << F2 << " result" << result <<std::endl;
	throw std::exception(); //("NaN in p_corr_R");
    }
    
    return result;
}
#endif

const Real 
BasicPairGreensFunction::p_corr_R( const Real alpha,
                                   const unsigned int n,
                                   const Real r,
                                   const Real r0,
                                   const Real t ) const
{
    const Real D( this->getD() );
    const Real sigma( this->getSigma() );
    
    const Real kSigma( getkf() * sigma );
    const Real realn( static_cast<Real>( n ) );
    const Real kSigma_m_n( kSigma - realn );

    const Real order( realn + 0.5 );

    const Real term1( exp( - D * t * alpha * alpha ) );

    const Real sigmaAlpha( sigma * alpha );
    const Real rAlpha( r * alpha );
    const Real r0Alpha( r0 * alpha );

    const Real js1( gsl_sf_bessel_Jnu( order,       sigmaAlpha ) );
    const Real ys1( gsl_sf_bessel_Ynu( order,       sigmaAlpha ) );
    const Real js2( gsl_sf_bessel_Jnu( order + 1.0, sigmaAlpha ) );
    const Real ys2( gsl_sf_bessel_Ynu( order + 1.0, sigmaAlpha ) );
    const Real jr(  gsl_sf_bessel_Jnu( order,       rAlpha ) );
    const Real yr(  gsl_sf_bessel_Ynu( order,       rAlpha ) );
    const Real jr0( gsl_sf_bessel_Jnu( order,       r0Alpha ) );
    const Real yr0( gsl_sf_bessel_Ynu( order,       r0Alpha ) );

    const Real R1( ( kSigma_m_n * js1 + sigmaAlpha * js2 ) * 2.0 );
    const Real R2( ( kSigma_m_n * ys1 + sigmaAlpha * ys2 ) * 2.0 );
    const Real F1R1( R1 * jr * jr0 - R1 * yr * yr0 );
    const Real F2( jr0 * yr + jr * yr0 );

    const Real num( alpha * R1 * ( F1R1 + F2 * R2 ) );
    const Real den( R1 * R1 + R2 * R2 );

    const Real result( term1 * num / den );

    //printf("%g %g %g %g %g %g\n", result, alpha, num, den, R1, R2 );

    assert( std::isfinite( result ) );

    return result;
}


const Real 
BasicPairGreensFunction::p_corr_R_F( const Real alpha, 
                                     const p_corr_R_params* const params )
{
    const BasicPairGreensFunction* const gf( params->gf ); 

    const unsigned int n( params->n );
    const Real r( params->r );
    const Real r0( params->r0 );
    const Real t( params->t );

    return gf->p_corr_R( alpha, n, r, r0, t );
}


const Real 
BasicPairGreensFunction::p_corr( const Real theta, const Real r, 
                                 const Real r0, const Real t ) const
{
    RealVector RnTable;
    makeRnTable( RnTable, r, r0, t );

    return p_corr_table( theta, r, r0, t, RnTable );
}

const Real 
BasicPairGreensFunction::ip_corr( const Real theta, const Real r, 
                                  const Real r0, const Real t ) const
{
    RealVector RnTable;
    makeRnTable( RnTable, r, r0, t );

    return ip_corr_table( theta, r, r0, t, RnTable );
}


const Real 
BasicPairGreensFunction::p_free( const Real theta, const Real r, const Real r0, 
				 const Real t ) const
{
    return p_theta_free( theta, r, r0, t, getD() );
}

const Real 
BasicPairGreensFunction::p_survival( const Real t, const Real r0 ) const
{
    return 1.0 - p_reaction( t, r0 );
}


const Real 
BasicPairGreensFunction::p_reaction( const Real t, const Real r0 ) const
{
    const Real kf( getkf() );
    const Real D( getD() );
    const Real sigma( getSigma() );
    const Real alpha( getalpha() );
    const Real kD( getkD() );

    return __p_reaction_irr( t, r0, kf, D, sigma, alpha, kD  );
}


const Real 
BasicPairGreensFunction::p_reaction_F( const Real t,
				       const p_reaction_params* const params )
{
    const BasicPairGreensFunction* const gf( params->gf ); 
    const Real kf( gf->getkf() );
    const Real D( gf->getD() );
    const Real sigma( gf->getSigma() );
    const Real alpha( gf->getalpha() );
    const Real kD( gf->getkD() );

    const Real r0( params->r0 );
    const Real rnd( params->rnd );

    return __p_reaction_irr( t, r0, kf, D, sigma, alpha, kD  ) - rnd;
}


const Real 
BasicPairGreensFunction::p_int_r( const Real r, 
                                  const Real t, 
                                  const Real r0 ) const
{
    const Real kf( getkf() );
    const Real D( getD() );
    const Real sigma( getSigma() );
    const Real alpha( getalpha() );
    const Real kD( getkD() );

    const Real Dt( D * t );

    const Real kf_kD( kf + kD );
    const Real Dt4( 4.0 * Dt );
    const Real sqrtDt4( sqrt( Dt4 ) );
    const Real ksigma2( 2.0 * kf * sigma );
    const Real alphasqrtt( alpha * sqrt( t ) );

    const Real r_r0__2s___sqrtDt4( ( r - 2.0 * sigma + r0 ) / sqrtDt4 );
    const Real r_r0__sqrtDt4( ( r - r0 ) / sqrtDt4 );
    const Real r0_s__sqrtDt4( ( r0 - sigma ) / sqrtDt4 );

    const Real term1( ( expm1( - gsl_pow_2( r_r0__2s___sqrtDt4  ) )
                        - expm1( - gsl_pow_2( r_r0__sqrtDt4 ) ) ) * 
                        sqrt( Dt / M_PI ) );

    const Real erf_r_r0__2s___sqrtDt4( erf( r_r0__2s___sqrtDt4 ) );
    const Real term2( kf_kD * r0 * erf( r_r0__sqrtDt4 ) 
                      + kf_kD * r0 * erf_r_r0__2s___sqrtDt4
                      + ksigma2 * 
                      ( erf( r0_s__sqrtDt4 ) - erf_r_r0__2s___sqrtDt4 ) );

    const Real term3( kf * sigma * W( r0_s__sqrtDt4, alphasqrtt ) 
                      - ( kf * r + kD * ( r - sigma ) ) *
                      W( r_r0__2s___sqrtDt4, alphasqrtt ) );

    const Real result( ( 1 / r0 ) * ( term1 + ( 1 / kf_kD ) * 
                                      ( ( 0.5 * term2 ) + term3 ) ) );

    return result;
}

const Real 
BasicPairGreensFunction::p_int_r_F( const Real r,
                                    const p_int_r_params* const params )
{
    const BasicPairGreensFunction* const gf( params->gf ); 

    const Real t( params->t );
    const Real r0( params->r0 );
    const Real rnd( params->rnd );

    return gf->p_int_r( r, t, r0 ) - rnd;
}


/*
const Real 
BasicPairGreensFunction::p_int_r_max( const Real t, 
                                      const Real r0 ) const
{
    const Real kf( getkf() );
    const Real D( getD() );
    const Real sigma( getSigma() );
    const Real alpha( getalpha() );
    const Real kD( getkD() );

    const Real Dt( D * t );

    const Real kf_kD( kf + kD );
    const Real Dt4( 4.0 * Dt );
    const Real sqrtDt4( sqrt( Dt4 ) );
    const Real alphasqrtt( alpha * sqrt( t ) );
    const Real kfsigma( kf * sigma );

    const Real r0_s__sqrtDt4( ( r0 - sigma ) / sqrtDt4 );

    const Real term1( kfsigma * erfc( r0_s__sqrtDt4 ) );
    const Real term2( kfsigma * W( r0_s__sqrtDt4, alphasqrtt ) );

    const Real den( kf_kD * r0 );

    const Real result( 1.0 - ( term1 + term2 ) / den );

    return result;
}
*/


const Real BasicPairGreensFunction::drawTime( const Real rnd, 
					      const Real r0 ) const
{
    const Real sigma( this->getSigma() );

    THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r0 >= sigma );

    Real low( 1e-100 );
    Real high( 100 );

    {
	const Real maxp( p_reaction( INFINITY, r0 ) );

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
  
    const Real r( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );

    return r;
} 




const Real BasicPairGreensFunction::drawR( const Real rnd, 
					   const Real r0, 
					   const Real t ) const
{
    const Real sigma( this->getSigma() );
    const Real D( this->getD() );

    THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r0 >= sigma );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    if( t == 0.0 )
    {
	return r0;
    }

    const Real psurv( p_survival( t, r0 ) );

    p_int_r_params params = { this, t, r0, rnd * psurv };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &p_int_r_F ),
	    &params 
	};


    // adjust low and high starting from r0.
    // this is necessary to avoid root finding in the long tails where
    // numerics can be unstable.

    Real low( r0 );
    Real high( r0 );

    const Real sqrt6Dt( sqrt( 6.0 * D * t ) );
    if( GSL_FN_EVAL( &F, r0 ) < 0.0 )
    {
        // low = r0
        unsigned int H( 3 );

        while( true )
        {
            high = r0 + H * sqrt6Dt;

            const Real value( GSL_FN_EVAL( &F, high ) );
            if( value > 0.0 )
            {
                break;
            }

            ++H;

            if( H > 20 )
            {
                std::cerr << "drawR: H > 20 while adjusting upper bound of r."
                          << std::endl;
                throw std::exception();
            }
        }

    }
    else
    {
        // high = r0
        unsigned int H( 3 );

        while( true )
        {
            low = r0 - H * sqrt6Dt;
            if( low < sigma )
            {
                if( GSL_FN_EVAL( &F, sigma ) > 0.0 )
                {
                    printf( "drawR: p_int_r( sigma ) > 0.0. "
                            "returning sigma.\n" );
                    return sigma;
                }

                low = sigma;
                break;
            }

            const Real value( GSL_FN_EVAL( &F, low ) );
            if( value < 0.0 )
            {
                break;
            }

            ++H;
        }
    }


    // root finding by iteration.

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
	const int status( gsl_root_test_interval( low, high, 1e-15,
						  this->TOLERANCE ) );

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


    const Real r( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );

    return r;
}


const Real 
BasicPairGreensFunction::Rn( const unsigned int n, const Real r, const Real r0,
			     const Real t,
			     gsl_integration_workspace* const workspace,
			     const Real tol ) const
{
    Real integral;
    Real error;

    p_corr_R_params params = { this, n, r, r0, t };
    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &p_corr_R_F ),
	    &params
	};

    const Real umax( sqrt( 40.0 / ( this->getD() * t ) ) ); 

    gsl_integration_qag( &F, 0.0,
			 umax,
			 tol,
			 1e-5,
			 2000, GSL_INTEG_GAUSS61,
			 workspace, &integral, &error );

/*
    gsl_integration_qagiu( &F, 1e-10,
                           0, // tol,
                           1e-6,
                           1000,
                           workspace, &integral, &error );
*/
  
    return integral;
}


const Real BasicPairGreensFunction::
p_corr_n( const unsigned int n, const RealVector& RnTable, 
          const RealVector& lgndTable ) const
{
    return RnTable[n] * lgndTable[n] * ( 2.0 * n + 1.0 );
}

const Real BasicPairGreensFunction::
ip_corr_n( const unsigned int n, const RealVector& RnTable, 
           const RealVector& lgndTable ) const
{
    // lgndTable1 is offset by 1; lgndTable1[0] is for n=-1.

    const Real lgnd_n_m1( lgndTable[n] );   // n-1
    const Real lgnd_n_p1( lgndTable[n+2] ); // n+1
    
    return RnTable[n] * ( lgnd_n_m1 - lgnd_n_p1 );// / ( 1.0 + 2.0 * n );
}


const Real BasicPairGreensFunction::
p_corr_table( const Real theta, const Real r, const Real r0,
              const Real t, const RealVector& RnTable ) const
{
    const Index tableSize( RnTable.size() );
    if( tableSize == 0 )
    {
        return 0.0;
    }

    Real result( 0.0 );


    Real sin_theta;
    Real cos_theta;
    sincos( theta, &sin_theta, &cos_theta );

    RealVector lgndTable( tableSize );
    gsl_sf_legendre_Pl_array( tableSize-1, cos( theta ), &lgndTable[0] );


    const Real p( funcSum_all( boost::bind( &BasicPairGreensFunction::
                                            p_corr_n,
                                            this,
                                            _1, RnTable, lgndTable ),
                               tableSize-1 ) );

    result = - p * sin_theta;

    result /= 4.0 * M_PI * sqrt( r * r0 );

    return result;
}


const Real BasicPairGreensFunction::
ip_corr_table( const Real theta, const Real r, const Real r0,
               const Real t, const RealVector& RnTable ) const
{
    const Index tableSize( RnTable.size() );
    if( tableSize == 0 )
    {
        return 0.0;
    }

    const Real cos_theta( cos( theta ) );
    
    // lgndTable is offset by 1. lengTable[0] -> n = -1

    RealVector lgndTable( tableSize + 2 );
    lgndTable[0] = 1.0; // n = -1
    gsl_sf_legendre_Pl_array( tableSize, cos_theta, &lgndTable[1] );

    const Real p( funcSum_all( boost::bind( &BasicPairGreensFunction::
                                            ip_corr_n,
                                            this,
                                            _1, RnTable, lgndTable ),
                               tableSize - 1 ) );

    const Real result( - p / ( 4.0 * M_PI * sqrt( r * r0 ) ) );
    return result;
}

const Real 
BasicPairGreensFunction::ip_free( const Real theta, 
                                  const Real r, 
                                  const Real r0, 
                                  const Real t ) const
{
    return ip_theta_free( theta, r, r0, t, getD() );
}


const Real BasicPairGreensFunction::
p_theta( const Real theta, const Real r, const Real r0, const Real t ) const
{
    RealVector RnTable;
    makeRnTable( RnTable, r, r0, t );

    return p_theta_table( theta, r, r0, t, RnTable );
}

const Real BasicPairGreensFunction::
ip_theta( const Real theta, const Real r, const Real r0, const Real t ) const
{
    RealVector RnTable;
    makeRnTable( RnTable, r, r0, t );

    return ip_theta_table( theta, r, r0, t, RnTable );
}


const Real BasicPairGreensFunction::
p_theta_table( const Real theta, const Real r, const Real r0,
               const Real t, const RealVector& RnTable ) const
{
    const Real p_free( this->p_free( theta, r, r0, t ) );
    const Real p_corr( this->p_corr_table( theta, r, r0, t, RnTable ) ); 

//    return p_free;
    return ( p_free + p_corr );
}

const Real BasicPairGreensFunction::
ip_theta_table( const Real theta, const Real r, const Real r0,
                const Real t, const RealVector& RnTable ) const
{
    const Real p_free( this->ip_free( theta, r, r0, t ) );
    const Real p_corr( this->ip_corr_table( theta, r, r0, t, RnTable ) ); 

    //printf("%g %g\n",p_free,p_corr);

//    return p_free;
    return ( p_free + p_corr );
}

static const Real p_free_max( const Real r,
                       const Real r0,
                       const Real t,
                       const Real D )
{
    const Real Dt4( 4.0 * D * t );
    const Real Dt4Pi( Dt4 * M_PI );

    const Real term1( exp( - gsl_pow_2( r - r0 ) / Dt4 ) );
    const Real term2( 1.0 / sqrt( Dt4Pi * Dt4Pi * Dt4Pi ) );

    return term1 * term2;
}

void BasicPairGreensFunction::makeRnTable( RealVector& RnTable,
                                           const Real r,
                                           const Real r0,
                                           const Real t ) const
{
    RnTable.clear();

    const Real sigma( getSigma() );
    const Real D( getD() );
    const Real kf( getkf() );

    {  
        // First, estimate the size of p_corr, and if it's small enough,
        // we don't need to calculate it in the first place.
        const Real pirr( p_irr( r, t, r0, kf, D, sigma ) );
        const Real ipfree_max( ip_free( M_PI, r, r0, t ) * 2 * M_PI * r * r );
        
        if( fabs( ( pirr - ipfree_max ) / ipfree_max ) < 1e-8 )
        {
            return;
        }
    }


    const Real pfreemax( p_free_max( r, r0, t, D ) );

    gsl_integration_workspace* 
        workspace( gsl_integration_workspace_alloc( 2000 ) );
    
    Real Rn_prev( 0.0 );
    const Real RnFactor( 1.0 / ( 4.0 * M_PI * sqrt( r * r0 ) ) );

    const Real integrationTolerance( pfreemax / RnFactor * 1e-6 );
    const Real truncationTolerance( pfreemax * 1e-7 );
    
    unsigned int n( 0 );
    while( true ) 
    {
        const Real Rn( this->Rn( n, r, r0, t, workspace, 
                                 integrationTolerance ) );
	
        RnTable.push_back( Rn );
	
        //printf("%d %g %g %g\n",n, Rn*RnFactor, pfreemax, integrationTolerance);
	
        // truncate when converged enough.
        if( fabs( Rn * RnFactor ) < fabs( truncationTolerance ) &&
            fabs( Rn ) < fabs( Rn_prev ) )
        {
            break;
        }

    
        if( n >= this->MAX_ORDER )
        {
            std::cerr << "Rn didn't converge." << std::endl;
            break;
        }
	
        Rn_prev = Rn;
	
        ++n;
    }
//    printf("%d \n",n);
    gsl_integration_workspace_free( workspace );
}

const Real BasicPairGreensFunction::ip_theta_F( const Real theta,
                                                const p_theta_params* params )
{
    const BasicPairGreensFunction* const gf( params->gf ); 
    const Real r( params->r );
    const Real r0( params->r0 );
    const Real t( params->t );
    const RealVector& RnTable( params->RnTable );
    const Real value( params->value );

    return gf->ip_theta_table( theta, r, r0, t, RnTable ) - value;
}


const Real BasicPairGreensFunction::drawTheta( const Real rnd,
					       const Real r, 
					       const Real r0, 
					       const Real t ) const
{
    Real theta;

    const Real sigma( this->getSigma() );

    // input parameter range checks.
    THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r >= sigma );
    THROW_UNLESS( std::invalid_argument, r0 >= sigma );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    // t == 0 means no move.
    if( t == 0.0 )
    {
	return 0.0;
    }

    RealVector RnTable;
    makeRnTable( RnTable, r, r0, t );


    // root finding with the integrand form.

    const Real ip_theta_pi( ip_theta_table( M_PI, r, r0, t, RnTable ) );

    p_theta_params params = { this, r, r0, t, RnTable, rnd * ip_theta_pi };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &ip_theta_F ),
	    &params 
	};

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, 0.0, M_PI );

    const unsigned int maxIter( 100 );

    unsigned int i( 0 );
    while( true )
    {
	gsl_root_fsolver_iterate( solver );
	const Real low( gsl_root_fsolver_x_lower( solver ) );
	const Real high( gsl_root_fsolver_x_upper( solver ) );
	const int status( gsl_root_test_interval( low, high, 1e-15,
						  this->TOLERANCE ) );

	if( status == GSL_CONTINUE )
	{
	    if( i >= maxIter )
	    {
		gsl_root_fsolver_free( solver );
		std::cerr << "drawTheta: failed to converge." << std::endl;
		throw std::exception();
	    }
	}
	else
	{
	    break;
	}

	++i;
    }
  
    theta = gsl_root_fsolver_root( solver );
    gsl_root_fsolver_free( solver );

    return theta;
}



//
// debug
//

const std::string BasicPairGreensFunction::dump() const
{
    std::ostringstream ss;
    ss << "D = " << this->getD() << ", sigma = " << this->getSigma() <<
	", kf = " << this->getkf() <<
	", kD = " << this->getkD() <<
	", alpha = " << this->getalpha() << std::endl;
    return ss.str();
}    

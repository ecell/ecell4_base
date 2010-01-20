#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <sstream>
#include <iostream>

#include <boost/format.hpp>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include "freeFunctions.hpp"

#include "FreePairGreensFunction.hpp"




const Real 
FreePairGreensFunction::p_r( const Real r, 
                             const Real r0, 
                             const Real t ) const
{
    const Real D( getD() );
    const Real Dt( D * t );
    const Real Dt4( 4.0 * Dt );
    const Real rr04( 4.0 * r * r0 );

    const Real mrr0sq_over_4Dt( - gsl_pow_2( r + r0 ) / Dt4 );

    const Real num1( expm1( mrr0sq_over_4Dt ) );
    const Real num2( expm1( mrr0sq_over_4Dt + rr04 / Dt4 ) );

    const Real den( rr04 * sqrt( M_PI * M_PI * M_PI * Dt ) );

    const Real jacobian( 2.0 * r * r * M_PI );

    return jacobian * ( - num1 + num2 ) / den;
}

const Real 
FreePairGreensFunction::ip_r( const Real r, 
                              const Real r0, 
                              const Real t ) const
{
    const Real D( getD() );
    const Real Dt4( 4.0 * D * t );
    const Real Dt4r( 1.0 / Dt4 );
    const Real sqrtDt4( sqrt( Dt4 ) );
    const Real sqrtDt4r( 1.0 / sqrtDt4 );

    const Real num1a( exp( - gsl_pow_2( r - r0 ) * Dt4r ) );
    const Real num1b( exp( - gsl_pow_2( r + r0 ) * Dt4r ) );
    const Real den1( r0 * sqrt( M_PI ) );

    const Real term1( sqrtDt4 * ( - num1a + num1b ) / den1 );

    const Real term2( erf( ( r - r0 ) * sqrtDt4r ) );
    const Real term3( erf( ( r + r0 ) * sqrtDt4r ) );

    return ( term1 + term2 + term3 ) * .5;
}
    
const Real 
FreePairGreensFunction::p_theta( const Real theta,
                                 const Real r,
                                 const Real r0,
                                 const Real t ) const
{
    return p_theta_free( theta, r, r0, t, getD() );
}

const Real 
FreePairGreensFunction::ip_theta( const Real theta,
                                  const Real r,
                                  const Real r0,
                                  const Real t ) const
{
    return ip_theta_free( theta, r, r0, t, getD() );
}

const Real
FreePairGreensFunction::ip_r_F( const Real r,
                                const ip_r_params* params )
{
    const FreePairGreensFunction* const gf( params->gf ); 
    const Real r0( params->r0 );
    const Real t( params->t );
    const Real value( params->value );

    return gf->ip_r( r, r0, t ) - value;
}



const Real 
FreePairGreensFunction::drawR( const Real rnd, 
                               const Real r0, 
                               const Real t ) const
{
    // input parameter range checks.
    if ( !(rnd <= 1.0 && rnd >= 0.0 ) )
    {
        throw std::invalid_argument( ( boost::format( "rnd <= 1.0 && rnd >= 0.0 : rnd=%g" ) % rnd ).str() );
    }

    if ( !(r0 >= 0.0 ) )
    {
        throw std::invalid_argument( ( boost::format( "r0 >= 0.0 : r0=%g" ) % r0 ).str() );
    }

    if ( !(t >= 0.0 ) )
    {
        throw std::invalid_argument( ( boost::format( "t >= 0.0 : t=%g" ) % t ).str() );
    }


    // t == 0 means no move.
    if( t == 0.0 )
    {
        return r0;
    }

    ip_r_params params = { this, r0, t, rnd };

    gsl_function F = 
        {
            reinterpret_cast<typeof(F.function)>( &ip_r_F ),
            &params 
        };

    const Real r_range( this->H * sqrt( 6.0 * getD() * t ) );

    const Real low_r( std::max( r0 - r_range, 0.0 ) );
    //const Real low_r( 0 );
    const Real max_r( r0 + r_range );


    if( GSL_FN_EVAL( &F, low_r ) >= 0.0 )
    {
        return low_r;
    }

    if( GSL_FN_EVAL( &F, max_r ) <= 0.0 )
    {
        return max_r;
    }

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, low_r, max_r );

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
FreePairGreensFunction::ip_theta_F( const Real theta,
                                    const ip_theta_params* params )
{
    const FreePairGreensFunction* const gf( params->gf ); 
    const Real r( params->r );
    const Real r0( params->r0 );
    const Real t( params->t );
    const Real value( params->value );

    return gf->ip_theta( theta, r, r0, t ) - value;
}



const Real 
FreePairGreensFunction::drawTheta( const Real rnd,
                                   const Real r, 
                                   const Real r0, 
                                   const Real t ) const
{
    // input parameter range checks.
    if ( !(rnd <= 1.0 && rnd >= 0.0 ) )
    {
        throw std::invalid_argument( ( boost::format( "rnd <= 1.0 && rnd >= 0.0 : rnd=%g" ) % rnd ).str() );
    }

    if ( !(r >= 0.0 ) )
    {
        throw std::invalid_argument( ( boost::format( "r >= 0.0 : r=%g" ) % r ).str() );
    }

    if ( !(r0 >= 0.0 ) )
    {
        throw std::invalid_argument( ( boost::format( "r0 >= 0.0 : r0=%g" ) % r0 ).str() );
    }

    if ( !(t >= 0.0 ) )
    {
        throw std::invalid_argument( ( boost::format( "t >= 0.0 : t=%g" ) % t ).str() );
    }


    // t == 0 means no move.
    if( t == 0.0 )
    {
        return 0.0;
    }

    const Real ip_theta_pi( ip_theta( M_PI, r, r0, t ) );

    ip_theta_params params = { this, r, r0, t, rnd * ip_theta_pi };

    gsl_function F = 
        {
            reinterpret_cast<typeof(F.function)>( &ip_theta_F ),
            &params 
        };

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, 0.0,
                          M_PI + std::numeric_limits<Real>::epsilon() );

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
  
    //printf("%d\n", i );

    const Real theta( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );
    

    return theta;
}


const std::string FreePairGreensFunction::dump() const
{
    std::ostringstream ss;
    ss << "D = " << this->getD() << std::endl;
    return ss.str();
}    

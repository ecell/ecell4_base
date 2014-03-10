#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <sstream>
#include <boost/format.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include "freeFunctions.hpp"
#include "GreensFunction3D.hpp"

GreensFunction3D::~GreensFunction3D()
{
    ; // do nothing
}
    
Real GreensFunction3D::drawTime(Real rnd) const
{
    return INFINITY;
}

Real GreensFunction3D::p_r(Real r, Real t) const
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

Real GreensFunction3D::ip_r(Real r, Real t) const
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
    
Real GreensFunction3D::p_theta(Real theta, Real r, Real t) const
{
    return p_theta_free( theta, r, r0, t, getD() );
}

Real GreensFunction3D::ip_theta(Real theta, Real r, Real t) const
{
    return ip_theta_free( theta, r, r0, t, getD() );
}

struct ip_r_params
{ 
    GreensFunction3D const* const gf;
    const Real t;
    const Real value;
};

static Real ip_r_F(Real r, ip_r_params const* params )
{
    return params->gf->ip_r(r, params->t) - params->value;
}

Real GreensFunction3D::drawR(Real rnd, Real t) const
{
    // input parameter range checks.
    if ( !(rnd <= 1.0 && rnd >= 0.0 ) )
    {
        throw std::invalid_argument( ( boost::format( "rnd <= 1.0 && rnd >= 0.0 : rnd=%.16g" ) % rnd ).str() );
    }

    if ( !(r0 >= 0.0 ) )
    {
        throw std::invalid_argument( ( boost::format( "r0 >= 0.0 : r0=%.16g" ) % r0 ).str() );
    }

    if ( !(t >= 0.0 ) )
    {
        throw std::invalid_argument( ( boost::format( "t >= 0.0 : t=%.16g" ) % t ).str() );
    }


    // t == 0 means no move.
    if( t == 0.0 )
    {
        return r0;
    }

    ip_r_params params = { this, t, rnd };

    gsl_function F = 
        {
            reinterpret_cast<typeof(F.function)>( &ip_r_F ),
            &params 
        };

    const Real r_range( this->H * sqrt( 6.0 * getD() * t ) );

    const Real low_r( std::max( r0 - r_range, 0.0 ) );
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
                throw std::runtime_error("drawR: failed to converge");
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

struct ip_theta_params
{ 
    GreensFunction3D const* const gf;
    const Real r;
    const Real t;
    const Real value;
};

static Real ip_theta_F(Real theta, ip_theta_params const* params)
{
    return params->gf->ip_theta(theta, params->r, params->t) - params->value;
}

Real GreensFunction3D::drawTheta(Real rnd, Real r, Real t) const
{
    // input parameter range checks.
    if ( !(rnd <= 1.0 && rnd >= 0.0 ) )
    {
        throw std::invalid_argument( ( boost::format( "rnd <= 1.0 && rnd >= 0.0 : rnd=%.16g" ) % rnd ).str() );
    }

    if ( !(r >= 0.0 ) )
    {
        throw std::invalid_argument( ( boost::format( "r >= 0.0 : r=%.16g" ) % r ).str() );
    }

    if ( !(r0 >= 0.0 ) )
    {
        throw std::invalid_argument( ( boost::format( "r0 >= 0.0 : r0=%.16g" ) % r0 ).str() );
    }

    if ( !(t >= 0.0 ) )
    {
        throw std::invalid_argument( ( boost::format( "t >= 0.0 : t=%.16g" ) % t ).str() );
    }


    // t == 0 means no move.
    if( t == 0.0 )
    {
        return 0.0;
    }

    const Real ip_theta_pi( ip_theta( M_PI, r, t ) );

    ip_theta_params params = { this, r, t, rnd * ip_theta_pi };

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
                throw std::runtime_error("drawTheta: failed to converge");
            }
        }
        else
        {
            break;
        }

        ++i;
    }
  
    const Real theta( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );
    

    return theta;
}


std::string GreensFunction3D::dump() const
{
    std::ostringstream ss;
    ss << "D = " << this->getD() << std::endl;
    return ss.str();
}

Logger& GreensFunction3D::log_(
    Logger::get_logger("GreensFunction3D"));


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdexcept>
#include <sstream>

#include <boost/format.hpp>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include "GreensFunction3DSym.hpp"

Real GreensFunction3DSym::p_r(Real r, Real t) const
{
    const Real D( getD() );
    const Real Dt( D * t );
    const Real Dt4( 4.0 * Dt );

    const Real Dt4Pi( Dt4 * M_PI );

    const Real term1( 1.0 / sqrt( gsl_pow_3( Dt4Pi ) ) );
    const Real term2( exp( - r * r / Dt4 ) );

    const Real jacobian( 4.0 * r * r * M_PI );

    return jacobian * term1 * term2;
}

Real GreensFunction3DSym::ip_r(Real r, Real t) const
{
    const Real D( getD() );
    const Real Dt( D * t );
    const Real sqrtDt_r( 1.0 / sqrt( D * t ) );
    const Real sqrtPi_r( 1.0 / sqrt( M_PI ) );

    const Real term1( exp( - r * r / ( 4.0 * Dt ) ) * 
                      r * sqrtDt_r * sqrtPi_r );
    const Real term2( erf( r * 0.5 * sqrtDt_r ) );

    return term2 - term1;
}

struct ip_r_params
{ 
    GreensFunction3DSym const* const gf;
    const Real t;
    const Real value;
};


static Real ip_r_F(Real r, ip_r_params const* params)
{
    const GreensFunction3DSym* const gf( params->gf ); 
    const Real t( params->t );
    const Real value( params->value );

    return gf->ip_r( r, t ) - value;
}


Real GreensFunction3DSym::drawR(Real rnd, Real t) const
{
    // input parameter range checks.
    if ( !(rnd <= 1.0 && rnd >= 0.0 ) )
    {
        throw std::invalid_argument( ( boost::format( "rnd <= 1.0 && rnd >= 0.0 : rnd=%.16g" ) % rnd ).str() );
    }

    if ( !(t >= 0.0 ) )
    {
        throw std::invalid_argument( ( boost::format( "t >= 0.0 : t=%.16g" ) % t ).str() );
    }


    // t == 0 or D == 0 means no move.
    if( t == 0.0 || getD() == 0.0 )
    {
        return 0.0;
    }

    ip_r_params params = { this, t, rnd };

    gsl_function F = 
        {
            reinterpret_cast<typeof(F.function)>( &ip_r_F ),
            &params 
        };

    Real max_r( 4.0 * sqrt( 6.0 * getD() * t ) );

    while( GSL_FN_EVAL( &F, max_r ) < 0.0 )
    {
        max_r *= 10;
    }

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, 0.0, max_r );

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
  
    //printf("%d\n", i );

    const Real r( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );
    
    return r;
}


std::string GreensFunction3DSym::dump() const
{
    std::ostringstream ss;
    ss << "D = " << this->getD() << std::endl;
    return ss.str();
}


Logger& GreensFunction3DSym::log_(
        Logger::get_logger("GreensFunction3DSym"));

//#define NDEBUG
//#define BOOST_DISABLE_ASSERTS

#include <iostream>
#include <stdexcept>
#include <vector>
#include <sstream>

#include <boost/bind.hpp>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
//#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sum.h>

#include "factorial.hpp"
#include "bessel.hpp"

//#include "HalfOrderBesselGenerator.hpp"

#include "FirstPassageNoCollisionPairGreensFunction.hpp"



FirstPassageNoCollisionPairGreensFunction::
FirstPassageNoCollisionPairGreensFunction( const Real D ) 
    :
    PairGreensFunction( D, 0, 0 ),
    a( INFINITY )
{
    ; // do nothing
}

FirstPassageNoCollisionPairGreensFunction::
~FirstPassageNoCollisionPairGreensFunction()
{
    ; // do nothing
}

void FirstPassageNoCollisionPairGreensFunction::seta( const Real a )
{
    this->a = a;

//    clearAlphaTable();
}


const Real 
FirstPassageNoCollisionPairGreensFunction::drawTime( const Real rnd, 
                                                     const Real r0 ) const
{

}

const Real 
FirstPassageNoCollisionPairGreensFunction::drawR( const Real rnd, 
                                                  const Real r0, 
                                                  const Real t ) const
{

}
    
const Real 
FirstPassageNoCollisionPairGreensFunction::drawTheta( const Real rnd,
                                                      const Real r, 
                                                      const Real r0, 
                                                      const Real t ) const
{

}

const Real p_survival_n( unsigned int n,
                         const Real t,
                         const Real r0 ) const
{
    const Real a( this->geta() );
    const Real D( this->getD() );

    const Real nPi_a( n * M_PI * _a );
        
    const Real term1( exp( - D * nPi_a * nPi_a * t ) );
    const Real term2( sin( nPi_a * r0 ) );
    const Real num( n * M_PI * r0 );
    
    Real value( a2 * term1 * term2 / num );
    if( n % 2 == 0 )
    {
        value = -value;
    }

    return value;
}

const Real p_survival( const Real t,
                       const Real r0 ) const
{
    Real result;

    const Real _a( 1.0 / a );
    const Real a2( 2.0 * a );
    const Dt ( this->getD() * t );

    const unsigned int MAXN( ceil( a / M_PI 
                                   * sqrt( log( this->TOLERANCE ) / Dt ) ) );

    const Real p1( p_survival_n( 1, t, r0 ) );
    const Real p1abs( fabs( p1 ) );
    result = p1;

    for( unsigned int n( 2 ); n < MAXN; ++n )
    {
        const Real p( p_survival_n( n, t, r0 ) );
        result += p;
    }


    return result;
}

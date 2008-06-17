#define BOOST_AUTO_TEST_MAIN

#define BOOST_TEST_MODULE SphericalBesselGenerator

#include <boost/mpl/list.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "SphericalBesselGenerator.hpp"

//#include "SphericalBesselGenerator.cpp"


const unsigned int maxn( 51 );
const unsigned int tableResolution( 300 );

const SphericalBesselGenerator& generator( getSphericalBesselGenerator() );


const Real TOLERANCE( 1e-5 );

BOOST_AUTO_TEST_CASE( testJ )
{
    const UnsignedInteger resolution( 100 );
    const Real maxz( std::max( 1000., static_cast<Real>( maxn * maxn ) ) * 2 );

    for( UnsignedInteger i( 0 ); i <= resolution; ++i )
    {
        const Real z( maxz * i / resolution );
        
        for( UnsignedInteger n( 0 ); n <= maxn; ++n )
        {
            const Real tj( generator.j( n, z ) );
            const Real j( gsl_sf_bessel_jl( n, z ) );
            
            BOOST_CHECK_CLOSE( j, tj, TOLERANCE );

            //printf("%d %g\n",n,z);
        }
    }

}


BOOST_AUTO_TEST_CASE( testY )
{
    const UnsignedInteger resolution( 100 );
    const Real maxz( std::max( 1000., static_cast<Real>( maxn * maxn ) ) * 2 );

    for( UnsignedInteger i( 1 ); i <= resolution; ++i )
    {
        const Real z( maxz * i / resolution );
        
        for( UnsignedInteger n( 0 ); n <= maxn; ++n )
        {
            const Real ty( generator.y( n, z ) );
            const Real y( gsl_sf_bessel_yl( n, z ) );
            
            BOOST_CHECK_CLOSE( y, ty, TOLERANCE );

            //printf("y %d %g\n",n,z);
        }
    }

}

/*
BOOST_AUTO_TEST_CASE( testJ_large_x )
{
    const UnsignedInteger resolution( 100 );
    const Real maxz( generator.maxz( maxn ) * 10000 );

    for( UnsignedInteger i( 0 ); i <= resolution; ++i )
    {
        const Real z( maxz * i / resolution );
        
        for( UnsignedInteger n( 0 ); n <= maxn; ++n )
        {
            const Real tj( generator.j( n, z ) );
            const Real j( gsl_sf_bessel_jl( n, z ) );
            
            BOOST_CHECK_CLOSE( j, tj, TOLERANCE );

            //printf("j large %d %g\n",n,z);
        }
    }

}
*/


#define BOOST_AUTO_TEST_MAIN

#define BOOST_TEST_MODULE SphericalBesselTable

#include <boost/mpl/list.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "SphericalBesselTable.hpp"


const unsigned int maxn( 50 );
const unsigned int tableResolution( 200 );

static SphericalBesselTable table( maxn, tableResolution );


BOOST_AUTO_TEST_CASE( testJ )
{
    for( unsigned int n( 0 ); n <= maxn; ++n )
    {
        const Real maxz( table.maxz( n ) * 5 );
        const unsigned int maxi( 500 );
        
        for( UnsignedInteger i( 0 ); i < maxi; ++i )
        {
            const Real z( i * maxz / maxi );

            const Real j( gsl_sf_bessel_jl( n, z ) );
            const Real tj( table.j( n, z ) );

            BOOST_CHECK_CLOSE( j, tj, 1e-7 );

            printf("%d\n",i);
        }
    }
        

}

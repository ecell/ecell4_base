#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#define BOOST_TEST_MODULE SphericalBesselTable

#define BOOST_AUTO_TEST_MAIN

#include <boost/mpl/list.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "SphericalBesselTable.hpp"


const unsigned int maxn( 50 );
const unsigned int tableResolution( 200 );

static SphericalBesselTable table( maxn, tableResolution );

const Real TOLERANCE( 1e-6 );


BOOST_AUTO_TEST_CASE( testJ )
{
    const UnsignedInteger resolution( 100 );
    const Real maxz( table.maxz( maxn ) * 1.1 );

    for( UnsignedInteger i( 0 ); i <= resolution; ++i )
    {
        const Real z( i * maxz / resolution );
        

        for( UnsignedInteger n( 0 ); n <= maxn; ++n )
        {
            const Real tj( table.j( n, z ) );
            const Real j( gsl_sf_bessel_jl( n, z ) );
            
            BOOST_CHECK_CLOSE( j, tj, TOLERANCE );
        }
    }

}

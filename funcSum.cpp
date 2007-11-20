#include <iostream>

#include <gsl/gsl_sum.h>

#include "funcSum.hpp"

const Real 
funcSum( boost::function<const Real( const unsigned int i )> f,
	 const size_t max_i,
	 const Real tolerance )
{
    Real sum( 0.0 );

    RealVector pTable;
    const Real p_0( f( 0 ) );
    if ( p_0 == 0.0 )
    {
	return 0.0;
    }

    const Real threshold( fabs( p_0 * tolerance ) );
    pTable.push_back( p_0 );

    bool extrapolationNeeded( true );

    RealVector::size_type i( 1 ); 
    while( i <= max_i )
    {
	const Real p_i( f( i ) );
	pTable.push_back( p_i );

	if( threshold >= fabs( p_i ) ) // '=' is important when p0 is so small.
	{
	    extrapolationNeeded = false;
	    break;
	}
	
	++i;
    }

    if( ! extrapolationNeeded )
    {
	sum = std::accumulate( pTable.begin(), pTable.end(), 0.0 );
    }
    else
    {
	// std::cerr << "Using series acceleration." << std::endl;
	Real error;
        /*
	gsl_sum_levin_u_workspace* 
	    workspace( gsl_sum_levin_u_alloc( i ) );
	gsl_sum_levin_u_accel( &pTable[0], pTable.size(), workspace, 
        &sum, &error );*/
	gsl_sum_levin_utrunc_workspace* 
	    workspace( gsl_sum_levin_utrunc_alloc( i ) );
	gsl_sum_levin_utrunc_accel( &pTable[0], pTable.size(), workspace, 
        &sum, &error );
	if( fabs( error ) >= fabs( sum * tolerance ) )
	{
	    std::cerr << "Series acceleration error; "
		      << fabs( error ) << " (rel error: " 
		      << fabs( error / sum ) << "), terms_used = " 
		      << workspace->terms_used << " (" 
		      << pTable.size() << " given)." << std::endl;
	}

//	gsl_sum_levin_u_free( workspace );
	gsl_sum_levin_utrunc_free( workspace );
    }

    return sum;
}

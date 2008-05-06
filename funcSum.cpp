#include <iostream>

#include <boost/bind.hpp>

#include <gsl/gsl_sum.h>

#include "funcSum.hpp"

const Real 
funcSum_all( boost::function<const Real( const unsigned int i )> f,
             const size_t max_i )
{
    Real sum( 0.0 );

    const Real p_0( f( 0 ) );
    if ( p_0 == 0.0 )
    {
	return 0.0;
    }

    sum = p_0;

    RealVector::size_type i( 1 ); 
    while( i <= max_i )
    {
	const Real p_i( f( i ) );
        sum += p_i;

	++i;
    }

    return sum;
}


const Real 
funcSum_all_accel( boost::function<const Real( const unsigned int i )> f,
                   const size_t max_i, const Real tolerance )
{
    RealVector pTable;
    pTable.reserve( max_i );

    const Real p_0( f( 0 ) );
    if ( p_0 == 0.0 )
    {
	return 0.0;
    }

    pTable.push_back( p_0 );

    RealVector::size_type i( 1 );
    for( ;  i <= max_i; ++i )
    {
	const Real p_i( f( i ) );
	pTable.push_back( p_i );
    }

    Real sum;
    Real error;
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
    
    gsl_sum_levin_utrunc_free( workspace );

    return sum;
}


const Real 
funcSum( boost::function<const Real( const unsigned int i )> f,
	 const size_t max_i,
	 const Real tolerance )
{
    const unsigned int CONVERGENCE_CHECK( 4 );

    Real sum( 0.0 );
    RealVector pTable;

    const Real p_0( f( 0 ) );
    if ( p_0 == 0.0 )
    {
	return 0.0;
    }

    //const Real threshold( fabs( p_0 * tolerance * 1e-1 ) );
    pTable.push_back( p_0 );
    sum = p_0;

    bool extrapolationNeeded( true );

    unsigned int convergenceCounter( 0 );

    RealVector::size_type i( 1 ); 
    while( i <= max_i )
    {
	const Real p_i( f( i ) );
	pTable.push_back( p_i );
        sum += p_i;

	++i;

	if( fabs( sum ) * tolerance >= fabs( p_i ) ) // '=' is important
        {
            ++convergenceCounter;
        }
        /*else  this screws it up; why?
        {
            convergenceCounter = 0;
            }*/

        if( convergenceCounter >= CONVERGENCE_CHECK )
	{
	    extrapolationNeeded = false;
	    break;
	}
	
    }

    if( extrapolationNeeded )
    {
        //std::cerr << "Using series acceleration." << i << std::endl;
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

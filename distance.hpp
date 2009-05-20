
#include <boost/multi_array.hpp>

#include <gsl/gsl_math.h>


const double _lengthSq( const double* const r )
{
    return gsl_pow_2( r[0] )
	+ gsl_pow_2( r[1] )
	+ gsl_pow_2( r[2] );
}

const double _length( const double* const r )
{
    return std::sqrt( _lengthSq( r ) );
}


const double _distanceSq( const double* const p1, const double* const p2 )
{
    return gsl_pow_2( p1[0] - p2[0] ) 
	+ gsl_pow_2( p1[1] - p2[1] ) 
	+ gsl_pow_2( p1[2] - p2[2] );
}

const double _distance( const double* const p1, const double* const p2 )
{
    return std::sqrt( _distanceSq( p1, p2 ) );
}


const double _distanceSq_Cyclic( const double* const p1, 
                                const double* const p2,
                                const double worldSize )
{
    const double halfWorldSize( worldSize * .5 );

    double diff[3] = { fabs( p1[0] - p2[0] ),
                       fabs( p1[1] - p2[1] ),
                       fabs( p1[2] - p2[2] ) };

    if( diff[0] > halfWorldSize )
    {
        diff[0] -= worldSize;
    }
    if( diff[1] > halfWorldSize )
    {
        diff[1] -= worldSize;
    }
    if( diff[2] > halfWorldSize )
    {
        diff[2] -= worldSize;
    }

    return gsl_pow_2( diff[0] )
        + gsl_pow_2( diff[1] )
	+ gsl_pow_2( diff[2] );
}

const double _distance_Cyclic( const double* const p1, 
                               const double* const p2,
                               const double worldSize )
{
    return std::sqrt( _distanceSq_Cyclic( p1, p2, worldSize ) );
}

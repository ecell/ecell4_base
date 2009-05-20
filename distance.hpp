#include <boost/multi_array.hpp>
#include <cmath>

double _lengthSq( double const (&r)[3] )
{
    return std::pow( r[0], 2) + std::pow( r[1], 2 ) + std::pow( r[2], 2 );
}

double _length( double const (&r)[3] )
{
    return std::sqrt( _lengthSq( r ) );
}


double _distanceSq( double const (&p1)[3], double const (&p2)[3] )
{
    return std::pow( p1[0] - p2[0], 2)
        + std::pow( p1[1] - p2[1], 2 ) 
        + std::pow( p1[2] - p2[2], 2 );
}

double _distance( double const (&p1)[3], double const (&p2)[3] )
{
    return std::sqrt( _distanceSq( p1, p2 ) );
}


double _distanceSq_Cyclic( double const (&p1)[3], 
                                double const (&p2)[3],
                                const double worldSize )
{
    const double halfWorldSize( worldSize * .5 );

    double diff[3] = { std::fabs( p1[0] - p2[0] ),
                       std::fabs( p1[1] - p2[1] ),
                       std::fabs( p1[2] - p2[2] ) };

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

    return std::pow( diff[0], 2) + std::pow( diff[1], 2 ) + std::pow( diff[2], 2 );
}

const double _distance_Cyclic( double const (&p1)[3], 
                               double const (&p2)[3],
                               const double worldSize )
{
    return std::sqrt( _distanceSq_Cyclic( p1, p2, worldSize ) );
}

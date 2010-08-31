#if !defined( __OLDDEFS_HPP )
#define __OLDDEFS_HPP

#include "config.h"
#include <vector>
#include <boost/multi_array.hpp>
#include <cmath>

// This file is needed temporarily by   GreensFunction1DAbsAbs.cpp
// and                                  GreensFunction1DRadAbs.cpp
//
// At some point it should be taken out of the new code version.
// This requires some major refurbishment of the above functions, so postponed for later.
//
// All passages conflicting with analogous definitions in the new Defs.hpp are taken out.
//




typedef std::vector< Real > RealVector;
typedef boost::multi_array<Real, 2>
Real2DArray;
typedef boost::multi_array<Real, 3>
Real3DArray;
typedef boost::multi_array<Real, 4>
Real4DArray;


#if !HAVE_SINCOS
inline void sincos( double x, double* s, double* c )
{
    *s = sin( x );
    *c = cos( x );
}
#endif /* !HAVE_SINCOS */

// stringifiers.  see preprocessor manual
#define XSTR( S ) STR( S )
#define STR( S ) #S

#define THROW_UNLESS( CLASS, EXPRESSION )       \
    if( ! ( EXPRESSION ) )\
    {\
        throw CLASS( "Check [" + std::string( STR( EXPRESSION ) ) +\
                     "] failed." );\
    }\


#define IGNORE_RETURN (void)


#endif // __OLDDEFS_HPP

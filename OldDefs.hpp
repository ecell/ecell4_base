#if !defined( __OLDDEFS_HPP )
#define __OLDDEFS_HPP

#include "config.h"
#include <vector>
#include <boost/multi_array.hpp>
#include <cmath>

// This file is needed temporarily by	FirstPassageGreensFunction1D.cpp
// and					FirstPassageGreensFunction1DRad.cpp
//
// At some point it should be taken out of the new code version.
// This requires some major refurbishment of the above functions, so postponed for later.
//
// All passages conflicting with analogous definitions in the new Defs.hpp are taken out.
//

// Taking out everything that is already defined in Defs.hpp
/*
//typedef double Real;
//typedef int Integer;
//typedef unsigned int UnsignedInteger;
typedef size_t Index;
*/


typedef std::vector< Real > RealVector;
//typedef boost::multi_array< Real, 1, boost::pool_allocator<Real> > 
//RealArray;
typedef boost::multi_array<Real, 2>
Real2DArray;
typedef boost::multi_array<Real, 3>
Real3DArray;
typedef boost::multi_array<Real, 4>
Real4DArray;

// Kick this out to prevent double declaration
/*
enum EventType
{
    REACTION = 0,
    ESCAPE = 1
};
*/


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

#define THROW_UNLESS( CLASS, EXPRESSION )	\
    if( ! ( EXPRESSION ) )\
    {\
	throw CLASS( "Check [" + std::string( STR( EXPRESSION ) ) +\
		     "] failed." );\
    }\


#define IGNORE_RETURN (void)


#endif // __OLDDEFS_HPP

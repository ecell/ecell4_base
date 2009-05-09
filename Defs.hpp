#if !defined( __DEFS_HPP )
#define __DEFS_HPP

#include "config.h"
#include <vector>
#include <boost/multi_array.hpp>
#include <cmath>

typedef double Real;
typedef long int Integer;
typedef unsigned long int UnsignedInteger;
typedef size_t Index;


typedef std::vector< Real > RealVector;
//typedef boost::multi_array< Real, 1, boost::pool_allocator<Real> > 
//RealArray;
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

#define THROW_UNLESS( CLASS, EXPRESSION )	\
    if( ! ( EXPRESSION ) )\
    {\
	throw CLASS( "Check [" + std::string( STR( EXPRESSION ) ) +\
		     "] failed." );\
    }\


#define IGNORE_RETURN (void)


#endif // __DEFS_HPP

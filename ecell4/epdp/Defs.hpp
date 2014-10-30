#if !defined( __DEFS_HPP )
#define __DEFS_HPP

#include <cstddef>
#include <ecell4/core/types.hpp>

// typedef double Real;
// typedef long int Integer;
// typedef unsigned long int UnsignedInteger;
typedef ecell4::Real Real;
typedef ecell4::Integer Integer;
typedef unsigned long int UnsignedInteger;
typedef size_t Index;

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

#endif // __DEFS_HPP

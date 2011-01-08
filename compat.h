#if !defined( COMPAT_HPP )
#define COMPAT_HPP

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#if defined( __cplusplus )
#include <cmath>
#include <limits>
#else
#include <math.h>
#include <limits.h>
#endif

#if !HAVE_DECL_INFINITY
#if defined( __cplusplus )
#    define INFINITY ( std::numeric_limits< double >::infinity() )
#else
#    if HAVE_DECL_HUGE_VAL
#        define INFINITY ( HUGE_VAL )
#    else
#        error could not define the constant 'INFINITY'
#    endif
#endif
#endif /* HAVE_DECL_INFINITY */

#if !defined( HAVE_SINCOS )
inline void sincos( double x, double* s, double* c )
{
    *s = sin( x );
    *c = cos( x );
}
#endif /* !HAVE_SINCOS */

#if !defined( HAVE_ISFINITE )
inline int isfinite( double x )
{
	return x == x && x != INFINITY && -x != INFINITY;
}
#endif

#endif // __COMPAT_HPP

#if !defined( __FUNCSUM_HPP )
#define __FUNCSUM_HPP

#include <boost/function.hpp>

#include "Defs.hpp"

static const Real TOLERANCE( 1e-6 );

const Real 
funcSum( boost::function<const Real( const unsigned int i )> f,
         const size_t max_i,
         const Real tolerance = TOLERANCE );


#endif /* __FUNCSUM_HPP */

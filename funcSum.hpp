#if !defined( __FUNCSUM_HPP )
#define __FUNCSUM_HPP

#include <boost/function.hpp>
#include <cstddef>

#include "Defs.hpp"

static const Real TOLERANCE( 1e-8 );

Real funcSum_all(boost::function<Real(unsigned int i)> f, std::size_t max_i);

Real funcSum_all_accel(boost::function<Real(unsigned int i)> f,
                       std::size_t max_i, Real tolerance = TOLERANCE);

Real funcSum(boost::function<Real(unsigned int i)> f,
             std::size_t max_i, Real tolerance = TOLERANCE);


#endif /* __FUNCSUM_HPP */

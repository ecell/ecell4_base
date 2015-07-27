#ifndef __ECELL4_TYPES_HPP
#define __ECELL4_TYPES_HPP

#include <stdint.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <limits>

namespace ecell4
{

typedef int64_t Integer;
typedef double Real;

const double inf = HUGE_VAL; // infinity (double)
const Real epsilon = std::numeric_limits<Real>::epsilon();
const Real N_A = 6.022140857e+23;

} // ecell4

#endif /* __ECELL4_TYPES_HPP */

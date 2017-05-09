#ifndef ECELL4_TYPES_HPP
#define ECELL4_TYPES_HPP

#include <stdint.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cfloat>
#include <limits>

namespace ecell4
{

typedef int64_t Integer;
typedef double Real;

const double inf = HUGE_VAL; // infinity (double)
const Real N_A = 6.022140857e+23;
const double epsilon = DBL_EPSILON; // std::numeric_limits<Real>::epsilon();

} // ecell4

#endif /* ECELL4_TYPES_HPP */

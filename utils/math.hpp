#ifndef UTILS_MATH_HPP
#define UTILS_MATH_HPP

#include <algorithm>
#include <cmath>

/**
 * Return True if a and b are equal, subject to given tolerances. Float 
 * comparison.
 *
 * See also numpy.allclose().
 *
 * The (relative) tolerance must be positive and << 1.0
 *
 * Instead of specifying an absolute tolerance, you can speciy a typical 
 * value for a or b. The absolute tolerance is then the relative tolerance 
 * multipied by this typical value, and will be used when comparing a value 
 * to zero. By default, the typical value is 1.
 */
template<typename T>
inline bool feq(T const& a, T const& b, T const& typical = 1., double tolerance = 1e-7)
{
    return std::abs(a - b) <= tolerance * (typical + std::min(std::abs(a), std::abs(b)));
}

#endif /* UTILS_MATH_HPP */

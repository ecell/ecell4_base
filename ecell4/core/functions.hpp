#ifndef __ECELL4_FUNCTIONS_HPP
#define __ECELL4_FUNCTIONS_HPP

#include <cmath>
#include <gsl/gsl_pow_int.h>

#include "config.h"
#include "types.hpp"


namespace ecell4
{

inline int64_t modulo(const int64_t& p1, const int64_t& p2)
{
    int64_t r(p1 % p2);
    if (r != 0 && (r > 0) == (p2 < 0))
    {
        r += p2;
    }
    return r;
}

inline double modulo(const double& p1, const double& p2)
{
    double r(std::fmod(p1, p2));
    if (r != 0 && (r > 0) == (p2 < 0))
    {
        r += p2;
    }
    return r;
}

inline double abs(const double& p1)
{
    return std::fabs(p1);
}

#ifndef WIN32_MSC
inline double pow_2(const double x)
{
    return gsl_pow_2(x);
}

inline double pow_3(const double x)
{
    return gsl_pow_3(x);
}
#else
inline double pow_2(const double x)
{
    return x * x;
}

inline double pow_3(const double x)
{
    return x * x * x;
}

inline double cbrt(const double x)
{
    return pow(x, 1.0 / 3.0);
}

inline double round(const double x)
{
    return floor(x + 0.5);
}
#endif

}

#endif /* __ECELL4_FUNCTIONS_HPP */

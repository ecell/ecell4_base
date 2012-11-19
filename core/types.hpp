#ifndef __TYPES_HPP
#define __TYPES_HPP

#include <stdint.h>
#include <cmath>


namespace ecell4
{

typedef int64_t Integer;
typedef double Real;

inline Integer modulo(Integer const& p1, Integer const& p2)
{
    Integer r = p1 % p2;
    if (r != 0 && (r > 0) == (p2 < 0))
        r += p2;
    return r;
}

inline Real modulo(Real const& p1, Real const& p2)
{
    Real r = std::fmod(p1, p2);
    if (r != 0 && (r > 0) == (p2 < 0))
        r += p2;
    return r;
}

inline Real abs(Real const& p1)
{
    return std::fabs(p1);
}

} // ecell4

#endif /* __TYPES_HPP */

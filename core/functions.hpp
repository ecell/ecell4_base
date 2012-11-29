#ifndef __FUNCTIONS_HPP
#define __FUNCTIONS_HPP

#include <cmath>


namespace ecell4
{

inline int64_t modulo(int64_t const& p1, int64_t const& p2)
{
    int64_t r(p1 % p2);
    if (r != 0 && (r > 0) == (p2 < 0))
    {
        r += p2;
    }
    return r;
}

inline double modulo(double const& p1, double const& p2)
{
    double r(std::fmod(p1, p2));
    if (r != 0 && (r > 0) == (p2 < 0))
    {
        r += p2;
    }
    return r;
}

inline double abs(double const& p1)
{
    return std::fabs(p1);
}

}

#endif /* __FUNCTIONS_HPP */

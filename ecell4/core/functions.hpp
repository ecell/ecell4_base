#ifndef __ECELL4_FUNCTIONS_HPP
#define __ECELL4_FUNCTIONS_HPP

#include <cmath>


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

}

#endif /* __ECELL4_FUNCTIONS_HPP */

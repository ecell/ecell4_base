#ifndef UTILS_RESET_HPP
#define UTILS_RESET_HPP

#include <boost/utility/swap.hpp>

template<typename T_>
inline void reset(T_& x)
{
    T_ y;
    boost::swap(y, x);
}

template<typename T_>
inline void reset(T_*& x)
{
    x = 0;
}

#endif /* UTILS_RESET_HPP */

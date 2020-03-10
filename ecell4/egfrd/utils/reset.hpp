#ifndef UTILS_RESET_HPP
#define UTILS_RESET_HPP

#include <algorithm>

namespace ecell4
{
namespace egfrd
{

template<typename T_>
inline void reset(T_& x)
{
    T_ y;
    std::swap(y, x);
}

template<typename T_>
inline void reset(T_*& x)
{
    x = 0;
}

} // egfrd
} // ecell4
#endif /* UTILS_RESET_HPP */

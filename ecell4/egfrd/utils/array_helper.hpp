#ifndef ECELL4_UTIL_ARRAY_HELPER_HPP
#define ECELL4_UTIL_ARRAY_HELPER_HPP
#include <array>

namespace ecell4
{
namespace egfrd
{

template<typename T>
inline std::array<T, 0> array_gen()
{
    return std::array<T, 0>{};
}
template<typename T>
inline std::array<T, 1> array_gen(const T& v0)
{
    return std::array<T, 1>{{v0}};
}
template<typename T>
inline std::array<T, 2> array_gen(const T& v0, const T& v1)
{
    return std::array<T, 2>{{v0, v1}};
}
template<typename T>
inline std::array<T, 3> array_gen(const T& v0, const T& v1, const T& v2)
{
    return std::array<T, 3>{{v0, v1, v2}};
}

} // egfrd
} // ecell4
#endif /* ARRAY_HELPER_HPP */

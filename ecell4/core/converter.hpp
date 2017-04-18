#ifndef ECELL4_UTIL_CONVERTER
#define ECELL4_UTIL_CONVERTER
#include <utility>

namespace ecell4
{
namespace utils
{

template<typename T_convert, typename T_first, typename T_second>
struct pair_first_element_converter
{
    typedef T_convert result_type;
    typedef std::pair<T_first, T_second> element_type;

    result_type operator()(const element_type& v) const
    {
        return result_type(v.first);
    }
};

template<typename T_convert, typename T_first, typename T_second>
struct pair_second_element_converter
{
    typedef T_convert result_type;
    typedef std::pair<T_first, T_second> element_type;

    result_type operator()(const element_type& v) const
    {
        return result_type(v.second);
    }
};

} // utils
} // ecell4
#endif// ECELL4_UTIL_CONVERTER

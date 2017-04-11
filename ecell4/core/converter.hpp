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
    typedef T_convert converted_type;
    typedef std::pair<T_first, T_second> element_type;

    converted_type operator()(const element_type& v)
    {
        return converted_type(v.first);
    }
};

template<typename T_convert, typename T_first, typename T_second>
struct pair_second_element_converter
{
    typedef T_convert converted_type;
    typedef std::pair<T_first, T_second> element_type;

    converted_type operator()(const element_type& v)
    {
        return converted_type(v.second);
    }
};

} // utils
} // ecell4
#endif// ECELL4_UTIL_CONVERTER

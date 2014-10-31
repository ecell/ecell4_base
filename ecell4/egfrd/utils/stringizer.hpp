#ifndef UTILS_STRINGIZER_HPP
#define UTILS_STRINGIZER_HPP

#include <string>
#include <functional>
#include <algorithm>
#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/value_type.hpp>
//#include "utils/range.hpp"
#include "./range.hpp"

template<typename T_>
struct stringizer: public std::unary_function<T_, std::string>
{
    std::string operator()(T_ const& value) const
    {
        return boost::lexical_cast<std::string>(value);
    }
};

template<typename T>
inline std::string stringize_and_join(T const& range, std::string const& separator)
{
    return boost::algorithm::join(
        make_transform_iterator_range(range,
            stringizer<typename boost::range_value<T>::type>()),
        separator);
}

#endif /* UTILS_STRINGIZER_HPP */

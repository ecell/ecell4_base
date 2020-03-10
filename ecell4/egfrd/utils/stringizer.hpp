#ifndef ECELL4_EGFRD_UTILS_STRINGIZER_HPP
#define ECELL4_EGFRD_UTILS_STRINGIZER_HPP

#include <string>
#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/adaptor/transformed.hpp>

namespace ecell4
{
namespace egfrd
{

template<typename T_>
struct stringizer
{
    std::string operator()(T_ const& value) const
    {
        return boost::lexical_cast<std::string>(value);
    }
};

template<typename T>
inline std::string stringize_and_join(T const& range, std::string const& separator)
{
    return boost::algorithm::join(range | boost::adaptors::transformed(
            stringizer<typename boost::range_value<T>::type>()), separator);
}

} // egfrd
} // ecell4
#endif /* ECELL4_EGFRD_UTILS_STRINGIZER_HPP */

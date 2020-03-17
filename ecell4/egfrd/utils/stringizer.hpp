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

template<typename Range>
std::string stringize_and_join(const Range& range, const std::string& separator)
{
    using value_type = typename boost::range_value<Range>::type;

    return boost::algorithm::join(boost::adaptors::transform(range,
                [](const value_type& v) -> std::string {
                    return boost::lexical_cast<std::string>(v);
                }), separator);
}

} // egfrd
} // ecell4
#endif /* ECELL4_EGFRD_UTILS_STRINGIZER_HPP */

#ifndef ABSTRACT_SET_HPP
#define ABSTRACT_SET_HPP

#include <iterator>
#include <algorithm>
#include <vector>
#include <set>
#include <map>

namespace ecell4
{
namespace egfrd
{

template<typename T_>
inline bool collection_contains(T_ const& s, typename T_::value_type const& v)
{
    auto e(std::end(s));
    return e != std::find(std::begin(s), e, v);
}

} // egfrd
} // ecell4
#endif /* ABSTRACT_SET_HPP */

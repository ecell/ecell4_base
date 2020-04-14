#ifndef ECELL4_EGFRD_COLLECTION_CONTAINS_HPP
#define ECELL4_EGFRD_COLLECTION_CONTAINS_HPP

#include <iterator>
#include <algorithm>

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
#endif /* ECELL4_EGFRD_COLLECTION_CONTAINS_HPP */

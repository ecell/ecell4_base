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

namespace detail
{

template<typename T, typename Alloc>
bool insert_impl(std::vector<T, Alloc>& set, const T& v)
{
    set.push_back(v);
    return true;
}

template<typename Key, typename Compare, typename Alloc>
bool insert_impl(std::set<Key, Compare, Alloc>& set, const Key& v)
{
    return set.insert(v).second;
}

template<typename Key, typename T, typename Compare, typename Alloc>
bool insert_impl(std::map<Key, T, Compare, Alloc>& set,
                 const typename std::map<Key, T, Compare, Alloc>::value_type& v)
{
    return set.insert(v).second;
}
} // detail

template<typename T_>
inline bool collection_contains(T_ const& s, typename T_::value_type const& v)
{
    auto e(std::end(s));
    return e != std::find(std::begin(s), e, v);
}

template<typename T_>
inline bool insert(T_& s, typename T_::value_type const& v)
{
    return detail::insert_impl(s, v);
}

template<typename T1, typename T2, typename OutputIterator>
inline void difference(T1 const& r1, T2 const& r2, OutputIterator result)
{
    std::set_difference(std::begin(r1), std::end(r1),
                        std::begin(r2), std::end(r2), result);
    return;
}

} // egfrd
} // ecell4
#endif /* ABSTRACT_SET_HPP */

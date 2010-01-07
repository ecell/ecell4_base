#ifndef RANGE_SUPPORT_HPP
#define RANGE_SUPPORT_HPP

#include <map>
#include <set>
#include <boost/range/size.hpp>
#include <boost/range/difference_type.hpp>

#if defined(HAVE_UNORDERED_MAP)
#include <unordered_map>
#endif

#if defined(HAVE_TR1_UNORDERED_MAP)
#include <tr1/unordered_map>
#endif

#if defined(HAVE_BOOST_UNORDERED_MAP_HPP)
#include <boost/unordered_map.hpp>
#endif

#define COMMA ,
#define SPECIALIZE_BOOST_SIZE(T) \
inline typename boost::range_difference<T>::type size(T const& r) \
{ \
    return r.size(); \
} \

namespace boost {

template<typename T1_, typename T2_, typename T3_, typename T4_>
SPECIALIZE_BOOST_SIZE(std::map<T1_ COMMA  T2_ COMMA  T3_ COMMA  T4_>)

template<typename T1_, typename T2_, typename T3_>
SPECIALIZE_BOOST_SIZE(std::set<T1_ COMMA  T2_ COMMA  T3_>)

#if defined(HAVE_UNORDERED_MAP)
template<typename T1_, typename T2_, typename T3_, typename T4_, typename T5_>
SPECIALIZE_BOOST_SIZE(std::unordered_map<T1_ COMMA  T2_ COMMA  T3_ COMMA  T4_ COMMA  T5_>)
#endif

#if defined(HAVE_TR1_UNORDERED_MAP)
template<typename T1_, typename T2_, typename T3_, typename T4_, typename T5_>
SPECIALIZE_BOOST_SIZE(std::tr1::unordered_map<T1_ COMMA  T2_ COMMA  T3_ COMMA  T4_ COMMA  T5_>)
#endif

#if defined(HAVE_BOOST_UNORDERED_MAP_HPP)
template<typename T1_, typename T2_, typename T3_, typename T4_, typename T5_>
SPECIALIZE_BOOST_SIZE(boost::unordered_map<T1_ COMMA  T2_ COMMA  T3_ COMMA  T4_ COMMA  T5_>)
#endif

} // namespace boost

#undef SPECIALIZE_BOOST_SIZE
#undef COMMA

#endif /* RANGE_SUPPORT_HPP */

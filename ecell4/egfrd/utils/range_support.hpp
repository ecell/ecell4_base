#ifndef RANGE_SUPPORT_HPP
#define RANGE_SUPPORT_HPP

#include <map>
#include <set>
#include <boost/range/size.hpp>
#include <boost/range/difference_type.hpp>
//#include "utils/range.hpp"
#include "./range.hpp"

#include <unordered_map>

namespace ecell4
{
namespace egfrd
{

#define COMMA ,
#define SPECIALIZE_BOOST_SIZE(T) \
inline typename boost::range_difference<T>::type size(T const& r) \
{ \
    return r.size(); \
} \

#define SPECIALIZE_RANGE_SIZE(T) \
struct range_size<T> \
{ \
    typedef typename T::size_type type; \
};

#define SPECIALIZE_RANGE_SIZE_RETRIEVER(T) \
struct range_size_retriever<T> \
{ \
    typedef typename range_size<T>::type result_type; \
    typedef T argument_type; \
\
    result_type operator()(argument_type const& range) const \
    { \
        return range.size(); \
    } \
};

} // egfrd
} // ecell4

namespace boost {

template<typename T1_, typename T2_, typename T3_, typename T4_>
SPECIALIZE_BOOST_SIZE(std::map<T1_ COMMA  T2_ COMMA  T3_ COMMA  T4_>)

template<typename T1_, typename T2_, typename T3_>
SPECIALIZE_BOOST_SIZE(std::set<T1_ COMMA  T2_ COMMA  T3_>)

template<typename T1_, typename T2_, typename T3_, typename T4_, typename T5_>
SPECIALIZE_BOOST_SIZE(std::unordered_map<T1_ COMMA  T2_ COMMA  T3_ COMMA  T4_ COMMA  T5_>)

} // namespace boost

namespace ecell4
{
namespace egfrd
{

template<typename T1_, typename T2_, typename T3_, typename T4_>
SPECIALIZE_RANGE_SIZE(std::map<T1_ COMMA  T2_ COMMA  T3_ COMMA  T4_>)

template<typename T1_, typename T2_, typename T3_, typename T4_>
SPECIALIZE_RANGE_SIZE_RETRIEVER(std::map<T1_ COMMA  T2_ COMMA  T3_ COMMA  T4_>)

template<typename T1_, typename T2_, typename T3_>
SPECIALIZE_RANGE_SIZE(std::set<T1_ COMMA  T2_ COMMA  T3_>)

template<typename T1_, typename T2_, typename T3_>
SPECIALIZE_RANGE_SIZE_RETRIEVER(std::set<T1_ COMMA  T2_ COMMA  T3_>)

template<typename T1_, typename T2_, typename T3_, typename T4_, typename T5_>
SPECIALIZE_RANGE_SIZE(std::unordered_map<T1_ COMMA  T2_ COMMA  T3_ COMMA  T4_ COMMA  T5_>)
template<typename T1_, typename T2_, typename T3_, typename T4_, typename T5_>
SPECIALIZE_RANGE_SIZE_RETRIEVER(std::unordered_map<T1_ COMMA  T2_ COMMA  T3_ COMMA  T4_ COMMA  T5_>)

#undef SPECIALIZE_BOOST_SIZE
#undef SPECIALIZE_RANGE_SIZE_RETRIEVER
#undef SPECIALIZE_RANGE_SIZE
#undef COMMA

} // egfrd
} // ecell4
#endif /* RANGE_SUPPORT_HPP */

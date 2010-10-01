#ifndef ABSTRACT_SET_HPP
#define ABSTRACT_SET_HPP

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/value_type.hpp>
#include <algorithm>

template<typename T_>
struct collection_value: public boost::range_value<T_>
{
};

template<typename T_>
inline bool contains(T_ const& s, typename collection_value<T_>::type const& v)
{
    typename boost::range_const_iterator<T_>::type e(boost::end(s));
    return e != std::find(boost::begin(s), e, v);
}

#endif /* ABSTRACT_SET_HPP */

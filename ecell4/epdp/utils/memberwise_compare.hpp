#ifndef MEMBERWISE_COMPARE_HPP
#define MEMBERWISE_COMPARE_HPP

#include <algorithm>
#include <utility>
#include <boost/range/size.hpp>
#include <boost/range/const_iterator.hpp>

template<typename Tlhs_, typename Trhs_>
inline int memberwise_compare(Tlhs_ const& lhs, Trhs_ const& rhs)
{
    typedef typename boost::range_const_iterator<Tlhs_>::type lhs_iterator;
    typedef typename boost::range_const_iterator<Trhs_>::type rhs_iterator;

    if (boost::size(lhs) <= boost::size(rhs))
    {
        std::pair<lhs_iterator, rhs_iterator> pair(
            std::mismatch(lhs.begin(), lhs.end(), rhs.begin()));
        if (pair.first == lhs.end())
            return boost::size(lhs) - boost::size(rhs);
        return *pair.first < *pair.second ?  -1:
                *pair.first > *pair.second ? 1: 0;
    }
    else if (boost::size(lhs) > boost::size(rhs))
    {
        std::pair<rhs_iterator, lhs_iterator> pair(
            std::mismatch(rhs.begin(), rhs.end(), lhs.begin()));
        if (pair.first == rhs.end())
            return 1;
        return *pair.first < *pair.second ? 1:
                *pair.first > *pair.second ? -1: 0;
    }
    return 0;
}

#endif /* MEMBERWISE_COMPARE_HPP */

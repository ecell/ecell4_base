#ifndef MEMBERWISE_COMPARE_HPP
#define MEMBERWISE_COMPARE_HPP

#include <algorithm>
#include <utility>
#include <boost/range/size.hpp>

template<typename Tlhs_, typename Trhs_>
inline int memberwise_compare(Tlhs_ const& lhs, Trhs_ const& rhs)
{
    if (boost::size(lhs) <= boost::size(rhs))
    {
        std::pair<typename Tlhs_::iterator, typename Trhs_::iterator> pair(
            std::mismatch(lhs.begin(), lhs.end(), rhs.begin()));
        if (pair.first == lhs.end())
            return boost::size(lhs) - boost::size(rhs);
        return *pair.first < *pair.second ?  -1:
                *pair.first > *pair.second ? 1: 0;
    }
    else if (boost::size(lhs) > boost::size(rhs))
    {
        std::pair<typename Trhs_::iterator, typename Tlhs_::iterator> pair(
            std::mismatch(rhs.begin(), rhs.end(), lhs.begin()));
        if (pair.first == rhs.end())
            return 1;
        return *pair.first < *pair.second ? 1:
                *pair.first > *pair.second ? -1: 0;
    }
    return 0;
}

#endif /* MEMBERWISE_COMPARE_HPP */

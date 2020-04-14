#ifndef ECELL4_EGFRD_UTILS_PAIR_HPP
#define ECELL4_EGFRD_UTILS_PAIR_HPP

#include <boost/range/value_type.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include "range.hpp"

// ---------------------------------------------------------------------------
// functors to extract first/second element of std::pair

namespace ecell4
{
namespace egfrd
{

template < typename T_ >
struct select_first
{
    typedef T_ argument_type;
    typedef typename T_::first_type result_type;

    typename T_::first_type& operator()( T_& pair ) const
    {
        return pair.first;
    }

    typename T_::first_type const& operator()( T_ const& pair ) const
    {
        return pair.first;
    }
};

template < typename T_ >
struct select_second
{
    typedef T_ argument_type;
    typedef typename T_::second_type result_type;

    typename T_::second_type& operator()( T_& pair ) const
    {
        return pair.second;
    }

    typename T_::second_type const& operator()( T_ const& pair ) const
    {
        return pair.second;
    }
};

// --------------------------------------------------------------------------
// helper aliases

template<typename Trange_>
using select_first_range_t = boost::transformed_range<
    select_first<typename boost::range_value<Trange_>::type>, Trange_>;

template<typename Trange_>
using select_second_range_t = boost::transformed_range<
    select_second<typename boost::range_value<Trange_>::type>, Trange_>;

// --------------------------------------------------------------------------
// make_select_first_range

template<typename Trange_>
inline select_first_range_t<Trange_>
make_select_first_range(Trange_& range)
{
    using value_type = typename boost::range_value<Trange_>::type;
    return boost::adaptors::transform(range, select_first<value_type>());
}

template<typename Trange_>
inline select_first_range_t<const Trange_>
make_select_first_range(Trange_ const& range)
{
    using value_type = typename boost::range_value<Trange_>::type;
    return boost::adaptors::transform(range, select_first<value_type>());
}

// --------------------------------------------------------------------------
// make_select_second_range

template<typename Trange_>
inline select_second_range_t<Trange_>
make_select_second_range(Trange_& range)
{
    using value_type = typename boost::range_value<Trange_>::type;
    return boost::adaptors::transform(range, select_second<value_type>());
}

template<typename Trange_>
inline select_second_range_t<const Trange_>
make_select_second_range(Trange_ const& range)
{
    using value_type = typename boost::range_value<Trange_>::type;
    return boost::adaptors::transform(range, select_second<value_type>());
}

} // egfrd
} // ecell4
#endif /* UTILS_PAIR_HPP */

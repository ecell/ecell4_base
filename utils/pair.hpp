#ifndef UTILS_PAIR_HPP
#define UTILS_PAIR_HPP

#include <boost/range/size.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/value_type.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/type_traits/remove_const.hpp>
#include "utils/range.hpp"

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

template<typename T_>
struct get_select_first_iterator
{
    typedef boost::transform_iterator<
        select_first<typename boost::iterator_value<T_>::type>, T_> type;
};

template<typename T_>
struct get_select_second_iterator
{
    typedef boost::transform_iterator<
        select_second<typename boost::iterator_value<T_>::type>, T_> type;
};

template<typename T_>
inline typename get_select_first_iterator<T_>::type
make_select_first_iterator(T_ const& iter)
{
    return typename get_select_first_iterator<T_>::type(iter,
        select_first<typename boost::iterator_value<T_>::type>());
        
}

template<typename T_>
inline typename get_select_second_iterator<T_>::type
make_select_second_iterator(T_ const& iter)
{
    return typename get_select_second_iterator<T_>::type(iter,
        select_second<typename boost::iterator_value<T_>::type>());
}

template<typename Trange_>
struct get_select_first_range
{
    typedef select_first<typename boost::range_value<Trange_>::type> functor_type;
    typedef typename get_transformed_range<Trange_, functor_type>::type type;
};

template<typename Trange_>
struct get_select_second_range
{
    typedef select_second<typename boost::range_value<Trange_>::type> functor_type;
    typedef typename get_transformed_range<Trange_, functor_type>::type type;
};

template<typename Trange_>
inline typename get_select_first_range<Trange_>::type
make_select_first_range(Trange_ const& range)
{
    typedef typename get_select_first_range<Trange_>::type type;
    return type(range, typename type::functor_type());
}

template<typename Trange_>
inline typename get_select_second_range<Trange_>::type
make_select_second_range(Trange_ const& range)
{
    typedef typename get_select_second_range<Trange_>::type type;
    return type(range, typename type::functor_type());
}

template<typename Tpair_>
struct remove_const_first
{
    typedef std::pair<typename boost::remove_const<typename Tpair_::first_type>::type,
                      typename Tpair_::second_type> type;
};

#endif /* UTILS_PAIR_HPP */

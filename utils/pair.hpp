#ifndef UTILS_PAIR_HPP
#define UTILS_PAIR_HPP
#include <boost/range/size.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/type_traits/remove_const.hpp>

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
        select_first<typename T_::value_type>, T_> type;
};

template<typename T_>
struct get_select_second_iterator
{
    typedef boost::transform_iterator<
        select_second<typename T_::value_type>, T_> type;
};

template<typename T_>
inline typename get_select_first_iterator<T_>::type
make_select_first_iterator(T_ const& iter)
{
    return typename get_select_first_iterator<T_>::type(iter,
        select_first<typename T_::value_type>());
        
}

template<typename T_>
inline typename get_select_second_iterator<T_>::type
make_select_second_iterator(T_ const& iter)
{
    return typename get_select_second_iterator<T_>::type(iter,
        select_second<typename T_::value_type>());
}

namespace detail {

template<typename Trange_>
struct select_first_range_impl
{
    typedef boost::iterator_range<typename get_select_first_iterator<
        typename boost::range_iterator<Trange_>::type >::type> type;
};

template<typename Trange_>
struct select_second_range_impl
{
    typedef boost::iterator_range<typename get_select_second_iterator<
        typename boost::range_iterator<Trange_>::type >::type> type;
};

template<typename Trange_>
struct select_first_range_impl<const Trange_>
{
    typedef boost::iterator_range<typename get_select_first_iterator<
        typename boost::range_const_iterator<Trange_>::type >::type> type;
};

template<typename Trange_>
struct select_second_range_impl<const Trange_>
{
    typedef boost::iterator_range<typename get_select_second_iterator<
        typename boost::range_const_iterator<Trange_>::type >::type> type;
};

} // namespace detail

template<typename Trange_>
class select_first_range: public detail::select_first_range_impl<Trange_>::type
{
    typedef typename detail::select_first_range_impl<Trange_>::type base_type;

public:
    select_first_range(Trange_ const& range)
        : base_type(boost::begin(range), boost::end(range)),
          size_(boost::size(range)) {}

    select_first_range(Trange_ const& range,
            typename base_type::size_type size)
        : base_type(boost::begin(range), boost::end(range)), size_(size) {}

    typename base_type::size_type size() const
    {
        return size_;
    }

private:
    typename base_type::size_type size_;
};

template<typename Trange_>
class select_second_range: public detail::select_second_range_impl<Trange_>::type
{
    typedef typename detail::select_second_range_impl<Trange_>::type base_type;

public:
    select_second_range(Trange_ const& range)
        : base_type(boost::begin(range), boost::end(range)),
          size_(boost::size(range)) {}

    select_second_range(Trange_ const& range,
            typename base_type::size_type size)
        : base_type(boost::begin(range), boost::end(range)), size_(size) {}

    typename base_type::size_type size() const
    {
        return size_;
    }

private:
    typename base_type::size_type size_;
};

template<typename Trange_>
inline select_first_range<Trange_>
make_select_first_range(Trange_ const& range)
{
    return select_first_range<Trange_>(range);
}

template<typename Trange_>
inline select_second_range<Trange_>
make_select_second_range(Trange_ const& range)
{
    return select_second_range<Trange_>(range);
}

template<typename Tpair_>
struct remove_const_first
{
    typedef std::pair<typename boost::remove_const<typename Tpair_::first_type>::type,
                      typename Tpair_::second_type> type;
};

namespace boost {

template<typename Trange_>
typename boost::range_difference<Trange_>::type size(select_first_range<Trange_> const& r)
{
    return r.size();
}

template<typename Trange_>
typename boost::range_difference<Trange_>::type size(select_second_range<Trange_> const& r)
{
    return r.size();
}

} // namespace boost

#endif /* UTILS_PAIR_HPP */

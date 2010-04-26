#ifndef UTILS_RANGE_HPP
#define UTILS_RANGE_HPP

#include <boost/range/size.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/iterator/iterator_categories.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/utility/enable_if.hpp>

template<typename Trange_, typename Tfun_>
struct get_transformed_range
{
    typedef typename boost::range_iterator<Trange_>::type iterator_type;
    typedef boost::transform_iterator<Tfun_, iterator_type> transform_iterator_type;
    typedef boost::iterator_range<transform_iterator_type> type;
};

template<typename Trange_, typename Tfun_>
struct get_transformed_range<const Trange_, Tfun_>
{
    typedef typename boost::range_const_iterator<Trange_>::type iterator_type;
    typedef boost::transform_iterator<Tfun_, iterator_type> transform_iterator_type;
    typedef boost::iterator_range<transform_iterator_type> type;
};


template<typename Trange_>
struct range_iterator_category
    : public boost::BOOST_ITERATOR_CATEGORY<typename boost::range_iterator<Trange_>::type>
{
};

template<typename Trange_, typename Ticat_>
struct check_range_iterator_category
    : public boost::is_convertible<
        typename boost::iterator_category_to_traversal<
            typename range_iterator_category<Trange_>::type >::type,
        typename boost::iterator_category_to_traversal<Ticat_>::type > {};

template<typename Trange_, typename Tfun_>
inline boost::iterator_range<
    boost::transform_iterator<
        Tfun_, typename boost::range_const_iterator<Trange_>::type> >
make_transform_iterator_range(Trange_ const& range, Tfun_ const& fun)
{
    typedef boost::transform_iterator<Tfun_,
        typename boost::range_const_iterator<Trange_>::type> transform_iterator;

    return boost::iterator_range<transform_iterator>(
        transform_iterator(boost::begin(range), fun),
        transform_iterator(boost::end(range), fun));
}

template<typename Trange_, typename Tfun_>
inline boost::iterator_range<
    boost::transform_iterator<
        Tfun_, typename boost::range_iterator<Trange_>::type> >
make_transform_iterator_range(Trange_& range, Tfun_ const& fun)
{
    typedef boost::transform_iterator<Tfun_,
        typename boost::range_iterator<Trange_>::type> transform_iterator;

    return boost::iterator_range<transform_iterator>(
        transform_iterator(boost::begin(range), fun),
        transform_iterator(boost::end(range), fun));
}

template<typename Titer_>
class sized_iterator_range: public boost::iterator_range<Titer_>
{
    typedef boost::iterator_range<Titer_> base_type;

public:
    sized_iterator_range(): size_(0) {}

    template<typename Taiter_>
    sized_iterator_range(Taiter_ begin, Taiter_ end)
        : base_type(begin, end), size_(end - begin) {}

    template<typename Trange_>
    sized_iterator_range(Trange_ const& r)
        : base_type(r), size_(boost::size(r)) {}

    template<typename Trange_>
    sized_iterator_range(Trange_& r)
        : base_type(r), size_(boost::size(r)) {}

    template<typename Taiter_>
    sized_iterator_range(Taiter_ begin, Taiter_ end,
            typename base_type::size_type size)
        : base_type(begin, end), size_(size) {}

    template<typename Trange_>
    sized_iterator_range(Trange_ const& r, typename base_type::size_type size)
        : base_type(r), size_(size) {}

    template<typename Trange_>
    sized_iterator_range(Trange_& r, typename base_type::size_type size)
        : base_type(r), size_(size) {}

    typename base_type::size_type size() const
    {
        return size_;
    }

private:
    typename base_type::size_type size_;
};

template<typename Tfn, typename Trange>
inline void call_with_size_if_randomly_accessible(
    Tfn& fn, Trange const &range,
    typename boost::enable_if<
        check_range_iterator_category<
            Trange, boost::random_access_traversal_tag> >::type* = 0)
{
    fn(boost::size(range));
}

template<typename Tfn, typename Trange>
inline void call_with_size_if_randomly_accessible(
    Tfn& fn, Trange const &range,
    typename boost::disable_if<
        check_range_iterator_category<
            Trange, boost::random_access_traversal_tag> >::type* = 0)
{
}

template<typename Tfn, typename Trange>
inline void call_with_size_if_randomly_accessible(
    Tfn const& fn, Trange const &range,
    typename boost::enable_if<
        check_range_iterator_category<
            Trange, boost::random_access_traversal_tag> >::type* = 0)
{
    fn(boost::size(range));
}

template<typename Tfn, typename Trange>
inline void call_with_size_if_randomly_accessible(
    Tfn const& fn, Trange const &range,
    typename boost::disable_if<
        check_range_iterator_category<
            Trange, boost::random_access_traversal_tag> >::type* = 0)
{
}

namespace boost {

template<typename Titer_>
typename boost::iterator_difference<Titer_>::type size(sized_iterator_range<Titer_> const& r)
{
    return r.size();
}

} // namespace boost

#endif /* UTILS_RANGE_HPP */

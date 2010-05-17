#ifndef UTILS_RANGE_HPP
#define UTILS_RANGE_HPP

#include <boost/range/size.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/not.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/iterator/iterator_categories.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/utility/enable_if.hpp>

template<typename Trange_>
struct range_iterator_category
    : boost::BOOST_ITERATOR_CATEGORY<typename boost::range_iterator<Trange_>::type> {};

template<typename Trange_, typename Ticat_>
struct check_range_iterator_category
    : boost::is_convertible<
        typename boost::iterator_category_to_traversal<
            typename range_iterator_category<Trange_>::type >::type,
        typename boost::iterator_category_to_traversal<Ticat_>::type > {};

template<typename Trange_>
struct is_sized
    : check_range_iterator_category<Trange_, boost::random_access_traversal_tag>
{
};

template<typename Trange_>
struct range_size: boost::range_difference<Trange_> {};

template<typename Trange_>
struct range_size_retriever
{
    typedef Trange_ argument_type;
    typedef typename range_size<Trange_>::type result_type;

    result_type operator()(argument_type const& range) const
    {
        return boost::size(range);
    }
};

template<typename Trange_>
inline typename range_size<Trange_>::type
size(Trange_ const& r)
{
    return range_size_retriever<Trange_>()(r);
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
        : base_type(r), size_(::size(r)) {}

    template<typename Trange_>
    sized_iterator_range(Trange_& r)
        : base_type(r), size_(::size(r)) {}

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

template<typename Trange_>
struct is_referencing_range: boost::mpl::false_ {};

template<typename Titer_>
struct is_referencing_range<std::pair<Titer_, Titer_> >: boost::mpl::true_ {};

template<typename Titer_>
struct is_referencing_range<boost::iterator_range<Titer_> >: boost::mpl::true_ {};

template<typename Titer_>
struct is_sized<sized_iterator_range<Titer_> >: boost::mpl::true_
{
};

template<typename Tfn, typename Trange>
inline void call_with_size_if_randomly_accessible(
    Tfn& fn, Trange const &range,
    typename boost::enable_if<is_sized<Trange> >::type* = 0)
{
    fn(::size(range));
}

template<typename Tfn, typename Trange>
inline void call_with_size_if_randomly_accessible(
    Tfn& fn, Trange const &range,
    typename boost::disable_if<is_sized<Trange> >::type* = 0)
{
}

template<typename Tfn, typename Trange>
inline void call_with_size_if_randomly_accessible(
    Tfn const& fn, Trange const &range,
    typename boost::enable_if<is_sized<Trange> >::type* = 0)
{
    fn(::size(range));
}

template<typename Tfn, typename Trange>
inline void call_with_size_if_randomly_accessible(
    Tfn const& fn, Trange const &range,
    typename boost::disable_if<is_sized<Trange> >::type* = 0)
{
}

template<typename Titer_>
struct range_size<sized_iterator_range<Titer_> >
{
    typedef std::size_t type;
};

template<typename Titer_>
struct range_size_retriever<sized_iterator_range<Titer_> >
{
    typedef sized_iterator_range<Titer_> argument_type;
    typedef typename range_size<argument_type>::type result_type;

    result_type operator()(argument_type const& range) const
    {
        return range.size();
    }
};

namespace detail {

template<typename Trange_ = void, bool N_ = boost::mpl::and_<
    boost::mpl::not_<boost::is_same<Trange_, void> >,
    is_sized<Trange_> >::value>
struct get_default_range_holder
{
    template<typename Titer_>
    struct apply
    {
        typedef boost::iterator_range<Titer_> type;
    };
};

template<typename Trange_>
struct get_default_range_holder<Trange_, true>
{
    template<typename Titer_>
    struct apply
    {
        typedef sized_iterator_range<Titer_> type;
    };
};

} // namespace detail

template<typename Titer_, typename Tfun_, typename Tholder_getter_ = detail::get_default_range_holder<> >
class transformed_range
{
public:
    typedef Titer_ base_iterator;
    typedef Tfun_ functor_type;
    typedef boost::transform_iterator<Tfun_, base_iterator> iterator;
    typedef typename Tholder_getter_::template apply<base_iterator>::type holder_type;
    typedef iterator const_iterator;
    typedef typename boost::iterator_value<base_iterator>::type value_type;
    typedef typename boost::iterator_difference<base_iterator>::type difference_type;
    typedef std::size_t size_type;
    typedef typename boost::iterator_reference<base_iterator>::type reference;
    typedef reference const_reference;

public:
    template<typename Trange>
    transformed_range(Trange const& range, Tfun_ const& fun)
        : base_(range), fun_(fun) {}

    holder_type base() const
    {
        return base_;
    }

    iterator begin() const
    {
        return iterator(boost::begin(base_), fun_);
    }

    iterator end() const
    {
        return iterator(boost::end(base_), fun_);
    }

    size_type size() const
    {
        return ::size(base_);
    }

    reference operator[](size_type idx) const
    {
        return fun_()(base_[idx]);
    }

private:
    holder_type base_;
    functor_type fun_;
};

template<typename Trange_, typename Tfun_>
struct get_transformed_range
{
    typedef transformed_range<typename boost::range_iterator<Trange_>::type, Tfun_, detail::get_default_range_holder<Trange_> > type;
};

template<typename Trange_, typename Tfun_>
struct get_transformed_range<const Trange_, Tfun_>
{
    typedef transformed_range<typename boost::range_const_iterator<Trange_>::type, Tfun_, detail::get_default_range_holder<Trange_> > type;
};

template<typename Trange_, typename Tfun_>
inline typename get_transformed_range<const Trange_, Tfun_>::type
make_transform_iterator_range(Trange_ const& range, Tfun_ const& fun)
{
    typedef typename get_transformed_range<const Trange_, Tfun_>::type transformed_range;
    return transformed_range(range, fun);
}

template<typename Trange_, typename Tfun_>
inline typename get_transformed_range<Trange_, Tfun_>::type
make_transform_iterator_range(Trange_ const& range, Tfun_ const& fun, bool)
{
    typedef typename get_transformed_range<Trange_, Tfun_>::type transformed_range;
    return transformed_range(range, fun);
}

template<typename Titer_, typename Tfun_, typename Tholder_getter_>
struct range_size<transformed_range<Titer_, Tfun_, Tholder_getter_> >: public range_size<typename Tholder_getter_::template apply<Titer_>::type> {};


template<typename Titer_, typename Tfun_, typename Tholder_getter_>
struct is_sized<transformed_range<Titer_, Tfun_, Tholder_getter_> >: public is_sized<typename Tholder_getter_::template apply<Titer_>::type> {};

template<typename Titer_, typename Tfun_, typename Tholder_getter_>
struct range_size_retriever<transformed_range<Titer_, Tfun_, Tholder_getter_> >
{
    typedef transformed_range<Titer_, Tfun_, Tholder_getter_> argument_type;
    typedef typename range_size<typename Tholder_getter_::template apply<Titer_>::type>::type result_type;

    result_type operator()(argument_type const& range) const
    {
        return range.size();
    }
};

#endif /* UTILS_RANGE_HPP */

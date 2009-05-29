#ifndef UTILS_HPP
#define UTILS_HPP

#include <map>
#include <vector>
#include <boost/range/size.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <algorithm>

namespace get_default_impl
{
    namespace std
    {
        template<typename Tkey_, typename Tval_>
        struct map
        {
            typedef ::std::map<Tkey_, Tval_> type;
        };

        template<typename Tval_>
        struct vector
        {
            typedef ::std::vector<Tval_> type; 
        };
    } // std
} // namespace get_default_impl

template < class T >
struct pointee_greater
{
    bool operator()( T x, T y ) const { return *y < *x; }
};


template < class T >
struct pointee_less
{
    bool operator()( T x, T y ) const { return *y > *x; }
};

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

namespace detail
{
    template < typename Tderived_, typename Tfun1_, typename Tfun2_,
               typename Tretval_ = typename Tfun1_::result_type >
    struct unary_compose_impl
    {
        typedef typename Tfun2_::argument_type argument_type;
        typedef typename Tfun1_::result_type result_type;

        unary_compose_impl( Tfun1_ const& f1, Tfun2_ const& f2 )
            : f1_( f1 ), f2_( f2 ) {}

        result_type operator()( argument_type const& val ) const
        {
            return f1_( f2_( val ) );
        }

        result_type operator()( argument_type const& val )
        {
            return f1_( f2_( val ) );
        }

        result_type operator()( argument_type& val ) const
        {
            return f1_( f2_( val ) );
        }

        result_type operator()( argument_type& val )
        {
            return f1_( f2_( val ) );
        }

    private:
        Tfun1_ f1_;
        Tfun2_ f2_;
    };

    template < typename Tderived_, typename Tfun1_, typename Tfun2_ >
    struct unary_compose_impl< Tderived_, Tfun1_, Tfun2_, void >
    {
        typedef typename Tfun2_::argument_type argument_type;
        typedef void result_type;

        unary_compose_impl( Tfun1_ const& f1, Tfun2_ const& f2 )
            : f1_( f1 ), f2_( f2 ) {}

        void operator()( argument_type const& val ) const
        {
            f1_( f2_( val ) );
        }

        void operator()( argument_type const& val )
        {
            f1_( f2_( val ) );
        }

        void operator()( argument_type& val ) const
        {
            f1_( f2_( val ) );
        }

        void operator()( argument_type& val )
        {
            f1_( f2_( val ) );
        }

    private:
        Tfun1_ f1_;
        Tfun2_ f2_;
    };
} // namespace detail

template < typename Tfun1_, typename Tfun2_ >
struct unary_compose
    : public detail::unary_compose_impl<unary_compose< Tfun1_, Tfun2_ >,
                                              Tfun1_, Tfun2_ >
{
public:
    unary_compose( Tfun1_ const& f1, Tfun2_ const& f2 )
        : detail::unary_compose_impl< unary_compose, Tfun1_, Tfun2_ >( f1, f2 ) {}
};

template < typename Tfun1_, typename Tfun2_ >
inline unary_compose< Tfun1_, Tfun2_ >
compose_unary( Tfun1_ const& f1, Tfun2_ const& f2 )
{
    return unary_compose< Tfun1_, Tfun2_ >( f1, f2 );
}

template < typename T_ >
struct delete_ptr
{
    typedef void result_type;
    typedef T_* argument_type;

    void operator()( T_* ptr )
    {
        delete ptr;
    }
};

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

template<typename Trange_>
inline boost::iterator_range<typename get_select_first_iterator<
    typename boost::range_iterator<Trange_>::type >::type>
make_select_first_range(Trange_ const& range)
{
    return boost::iterator_range<typename get_select_first_iterator<
        typename boost::range_iterator<Trange_>::type >::type>(
            boost::begin(range), boost::end(range));
}

template<typename Trange_>
inline boost::iterator_range<typename get_select_second_iterator<
    typename boost::range_iterator<Trange_>::type >::type>
make_select_second_range(Trange_ const& range)
{
    return boost::iterator_range<typename get_select_second_iterator<
        typename boost::range_iterator<Trange_>::type >::type>(
            boost::begin(range), boost::end(range));
}

void gsl_error_handler( char const* reason, char const* file, int line, int gsl_errno );

#endif /* UTILS_HPP */

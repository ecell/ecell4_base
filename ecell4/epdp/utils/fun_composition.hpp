#ifndef FUN_COMPOSITION_HPP
#define FUN_COMPOSITION_HPP

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

    template < typename Tderived_, typename Tfun1_, typename Tfun2_,
               typename Tretval_ = typename Tfun1_::result_type >
    struct binary_compose_impl
    {
        typedef typename Tfun2_::first_argument_type first_argument_type;
        typedef typename Tfun2_::second_argument_type second_argument_type;
        typedef typename Tfun1_::result_type result_type;

        binary_compose_impl( Tfun1_ const& f1, Tfun2_ const& f2 )
            : f1_( f1 ), f2_( f2 ) {}

        result_type operator()( first_argument_type const& v1,
                                second_argument_type const& v2 ) const
        {
            return f1_( f2_( v1, v2 ) );
        }

        result_type operator()( first_argument_type const& v1,
                                second_argument_type const& v2 )
        {
            return f1_( f2_( v1, v2 ) );
        }

        result_type operator()( first_argument_type& v1,
                                second_argument_type& v2 ) const
        {
            return f1_( f2_( v1, v2 ) );
        }

        result_type operator()( first_argument_type& v1,
                                second_argument_type& v2 )
        {
            return f1_( f2_( v1, v2 ) );
        }

    private:
        Tfun1_ f1_;
        Tfun2_ f2_;
    };

    template < typename Tderived_, typename Tfun1_, typename Tfun2_ >
    struct binary_compose_impl< Tderived_, Tfun1_, Tfun2_, void >
    {
        typedef typename Tfun2_::first_argument_type first_argument_type;
        typedef typename Tfun2_::second_argument_type second_argument_type;
        typedef void result_type;

        binary_compose_impl( Tfun1_ const& f1, Tfun2_ const& f2 )
            : f1_( f1 ), f2_( f2 ) {}

        void operator()( first_argument_type const& v1,
                                second_argument_type const& v2 ) const
        {
            f1_( f2_( v1, v2 ) );
        }

        void operator()( first_argument_type const& v1,
                                second_argument_type const& v2 )
        {
            f1_( f2_( v1, v2 ) );
        }

        void operator()( first_argument_type& v1,
                                second_argument_type& v2 ) const
        {
            f1_( f2_( v1, v2 ) );
        }

        void operator()( first_argument_type& v1,
                                second_argument_type& v2 )
        {
            f1_( f2_( v1, v2 ) );
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

template < typename Tfun1_, typename Tfun2_ >
struct binary_compose
    : public detail::binary_compose_impl<binary_compose< Tfun1_, Tfun2_ >,
                                              Tfun1_, Tfun2_ >
{
public:
    binary_compose( Tfun1_ const& f1, Tfun2_ const& f2 )
        : detail::binary_compose_impl< binary_compose, Tfun1_, Tfun2_ >( f1, f2 ) {}
};

template < typename Tfun1_, typename Tfun2_ >
inline binary_compose< Tfun1_, Tfun2_ >
compose_binary( Tfun1_ const& f1, Tfun2_ const& f2 )
{
    return binary_compose< Tfun1_, Tfun2_ >( f1, f2 );
}

#endif /* FUN_COMPOSITION_HPP */

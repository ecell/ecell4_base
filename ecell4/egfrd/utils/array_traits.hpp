#ifndef ARRAY_CAST_HPP
#define ARRAY_CAST_HPP

#include <boost/array.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/mpl/if.hpp>
#include <boost/multi_array.hpp>

template< typename T_ >
struct num_elements
{
    struct cannot_deduce_number_of_elements_from_the_specified_type;
    enum { _ = sizeof(cannot_deduce_number_of_elements_from_the_specified_type) };
};

template< typename T_, std::size_t N_ >
struct num_elements< T_[ N_ ] >
{
    BOOST_STATIC_CONSTANT( std::size_t, value = N_ );
};

template< typename T_ >
struct element_type_of
{
//     typedef typename T_::value_type type;
};

template< typename T_, std::size_t N_ >
struct element_type_of< T_[N_] >
{
    typedef T_ type;
};

template< typename T_, std::size_t N_ >
struct element_type_of< boost::array< T_, N_ > >
{
    typedef T_ type;
};

template< typename T_, typename Talloc_ >
struct element_type_of< boost::multi_array< T_, 1, Talloc_ > >
{
    typedef T_ type;
};

template< typename T_ >
T_& array_cast( typename element_type_of< T_ >::type* elts )
{
    return *elts;
}

template< typename T_, typename Tx_ >
T_& array_cast( Tx_& elts )
{
    BOOST_STATIC_ASSERT(( boost::is_same<
        typename boost::remove_const<
            typename element_type_of< T_ >::type >::type,
        typename Tx_::value_type >::value ));
    return reinterpret_cast< T_& >( *elts.data() );
}

#endif /* ARRAY_CAST_HPP */

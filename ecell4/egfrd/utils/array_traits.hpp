#ifndef ARRAY_CAST_HPP
#define ARRAY_CAST_HPP

#include <boost/array.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/mpl/if.hpp>
#include <boost/multi_array.hpp>
namespace ecell4
{
namespace egfrd
{
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
struct element_type_of< std::array< T_, N_ > >
{
    typedef T_ type;
};

template< typename T_, typename Talloc_ >
struct element_type_of< boost::multi_array< T_, 1, Talloc_ > >
{
    typedef T_ type;
};

} // egfrd
} // ecell4
#endif /* ARRAY_CAST_HPP */

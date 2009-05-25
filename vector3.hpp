#ifndef VECTOR3_HPP
#define VECTOR3_HPP

#include <ostream>
#include <functional>
#include <algorithm>
#include <cmath>
#include <boost/array.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>

#include "array_traits.hpp"

template< typename T_ >
inline T_ add( T_ const& p1, T_ const& p2 )
{
    T_ retval;
    retval[0] = p1[0] + p2[0];
    retval[1] = p1[1] + p2[1];
    retval[2] = p1[2] + p2[2];
    return retval;
}

template< typename T_ >
inline T_ subtract( T_ const& p1, T_ const& p2 )
{
    T_ retval;
    retval[0] = p1[0] - p2[0];
    retval[1] = p1[1] - p2[1];
    retval[2] = p1[2] - p2[2];
    return retval;
}

template< typename T_ >
inline typename detail::element_type_of< T_ >::type length_sq( T_ const& r )
{
    return std::pow( r[0], 2) + std::pow( r[1], 2 ) + std::pow( r[2], 2 );
}

template< typename T_ >
inline typename detail::element_type_of< T_ >::type length( T_ const& r )
{
    return std::sqrt( length_sq( r ) );
}

template< typename T1_, typename T2_ >
inline typename detail::element_type_of< T1_ >::type distance_sq(
        T1_ const& p1, T2_ const p2 )
{
    return std::pow( p1[0] - p2[0], 2)
        + std::pow( p1[1] - p2[1], 2 ) 
        + std::pow( p1[2] - p2[2], 2 );
}

template< typename T1_, typename T2_ >
inline typename detail::element_type_of< T1_ >::type distance(
        T1_ const& p1, T2_ const& p2 )
{
    return std::sqrt( distance_sq( p1, p2 ) );
}

template< typename T1_, typename T2_ >
inline typename detail::element_type_of< T1_ >::type
distance_sq_cyclic( T1_ const& p1, T2_ const& p2,
        typename element_type_of< T1_ >::type const& world_size )
{
    typedef typename element_type_of< T1_ >::type element_type;
    const element_type half_world_size( world_size * .5 );

    element_type diff[3] = {
        std::fabs( p1[0] - p2[0] ),
        std::fabs( p1[1] - p2[1] ),
        std::fabs( p1[2] - p2[2] )
    };

    if( diff[0] > half_world_size )
    {
        diff[0] -= world_size;
    }

    if( diff[1] > half_world_size )
    {
        diff[1] -= world_size;
    }

    if( diff[2] > half_world_size )
    {
        diff[2] -= world_size;
    }

    return std::pow( diff[0], 2 ) +
        std::pow( diff[1], 2 ) +
        std::pow( diff[2], 2 );
}

template< typename T1_,typename T2_ >
typename detail::element_type_of< T1_ >::type
distance_cyclic( T1_ const& p1, T2_ const& p2,
        typename detail::element_type_of< T1_ >::type const& world_size )
{
    return std::sqrt( distance_sq_cyclic( p1, p2, world_size ) );
}

template< typename T1_ >
inline typename detail::element_type_of< T1_ >::type
distance( T1_ const& p1, T1_ const& p2 )
{
    return distance< T1_, T1_ >( p1, p2 );
}

template< typename T1_ >
inline typename detail::element_type_of< T1_ >::type
distance_sq( T1_ const& p1, T1_ const& p2 )
{
    return distance_sq< T1_, T1_ >( p1, p2 );
}

template< typename T1_ >
inline typename detail::element_type_of< T1_ >::type
distance_sq_cyclic( T1_ const& p1, T1_ const& p2,
        typename element_type_of< T1_ >::type const& world_size )
{
    return distance_sq_cyclic< T1_, T1_ >( p1, p2, world_size );
}

template< typename T1_ >
inline typename detail::element_type_of< T1_ >::type
distance_cyclic( T1_ const& p1, T1_ const& p2,
        typename element_type_of< T1_ >::type const& world_size )
{
    return distance_cyclic< T1_, T1_ >( p1, p2, world_size );
}


template<typename T_>
struct vector3: public boost::array<T_, 3>
{
    typedef boost::array<T_, 3> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::size_type size_type;

    vector3()
    {
        (*this)[0] = 0;
        (*this)[1] = 0;
        (*this)[2] = 0;
    }

    vector3(const T_ (&a)[3]): base_type(
            *reinterpret_cast<const base_type*>(&a)) {}

    vector3(const T_ a[3]): base_type(
            *reinterpret_cast<const base_type*>(a)) {}

    vector3(const base_type& a): base_type(a) {}

    vector3(value_type p0, value_type p1, value_type p2)
    {
        (*this)[0] = p0;
        (*this)[1] = p1;
        (*this)[2] = p2;
    }
};

template< typename T_ >
vector3< T_ > operator+(vector3< T_ > const& lhs, vector3< T_ > const& rhs)
{
    return add( lhs, rhs );
}

template< typename T_ >
vector3< T_ > operator-(vector3< T_ > const& lhs, vector3< T_ > const& rhs)
{
    return subtract( lhs, rhs );
}

template<typename Tstrm_, typename T_>
inline std::basic_ostream<Tstrm_>&
operator<<(std::basic_ostream< Tstrm_ >& strm, const vector3< T_ >& v)
{
    strm << "(" << v[0] <<  ", " << v[1] <<  ", " << v[2] << ")";
    return strm;
}

namespace detail {
    template< typename T_ >
    struct element_type_of< vector3< T_ > >
    {
        typedef T_ type;
    };
}

#endif /* VECTOR3_HPP */

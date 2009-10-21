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

template< typename T1_ >
inline T1_ normalize( T1_ const& p, 
                      typename detail::element_type_of< T1_ >::type const& r )
{
    typedef typename element_type_of< T1_ >::type element_type;

    T1_ retval;
    element_type factor(r / length<T1_>(p));
    retval[0] = p[0] * factor;
    retval[1] = p[1] * factor;
    retval[2] = p[2] * factor;

    return retval;
}

// cyclic_transpose() and calculatePairCoM are placed here
// temporarily.  should be moved to somewhere else later.

template< typename T1_ >
inline T1_ 
cyclic_transpose( T1_ const& p1, 
                  T1_ const& p2, 
                  typename detail::element_type_of< T1_ >::type const& 
                  world_size )
{
    T1_ retval;

    typedef typename element_type_of< T1_ >::type element_type;   
    const element_type half_world_size(world_size * .5);

    for(unsigned int i(0); i <= 2; ++i)
    {
        const element_type diff(p2[i] - p1[i]);
        retval[i] = p1[i];
        if(diff > half_world_size)
        {
            retval[i] += world_size;
        }
        else if(diff < - half_world_size)
        {
            retval[i] -= world_size;
        }
    }

    return retval;
}

template< typename T1_ >
inline T1_ 
calculate_pair_CoM( T1_ const& p1, 
                    T1_ const& p2, 
                    typename detail::element_type_of< T1_ >::type const& D1,
                    typename detail::element_type_of< T1_ >::type const& D2,
                    typename detail::element_type_of< T1_ >::type const& 
                    world_size )
{
    typedef typename element_type_of< T1_ >::type element_type;   

    T1_ retval;

    const T1_ p2t(cyclic_transpose<T1_>(p2,p1,world_size));

    const element_type rD1pD2(1.0/(D1 + D2));
    const element_type fD1(D1 * rD1pD2);
    const element_type fD2(D2 * rD1pD2);
    for(unsigned int i(0); i <= 2; ++i)
    {
        retval[i] = fmod((fD2 * p1[i] + fD1 * p2t[i]), world_size);
    }
    
    return retval;
}




template<typename T_>
struct Vector3: public boost::array<T_, 3>
{
    typedef boost::array<T_, 3> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::size_type size_type;

    Vector3()
    {
        (*this)[0] = 0;
        (*this)[1] = 0;
        (*this)[2] = 0;
    }

    Vector3(const T_ (&a)[3]): base_type(
            *reinterpret_cast<const base_type*>(&a)) {}

    Vector3(const T_ a[3]): base_type(
            *reinterpret_cast<const base_type*>(a)) {}

    Vector3(const base_type& a): base_type(a) {}

    Vector3(value_type p0, value_type p1, value_type p2)
    {
        (*this)[0] = p0;
        (*this)[1] = p1;
        (*this)[2] = p2;
    }
};

template< typename T_ >
Vector3< T_ > operator+(Vector3< T_ > const& lhs, Vector3< T_ > const& rhs)
{
    return add( lhs, rhs );
}

template< typename T_ >
Vector3< T_ > operator-(Vector3< T_ > const& lhs, Vector3< T_ > const& rhs)
{
    return subtract( lhs, rhs );
}

template<typename Tstrm_, typename Ttraits_, typename T_>
inline std::basic_ostream<Tstrm_, Ttraits_>&
operator<<(std::basic_ostream<Tstrm_, Ttraits_>& strm, const Vector3<T_>& v)
{
    strm << "(" << v[0] <<  ", " << v[1] <<  ", " << v[2] << ")";
    return strm;
}

namespace detail {
    template< typename T_ >
    struct element_type_of< Vector3< T_ > >
    {
        typedef T_ type;
    };
}

#endif /* VECTOR3_HPP */

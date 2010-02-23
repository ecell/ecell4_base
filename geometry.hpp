#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/and.hpp>
#include "linear_algebra.hpp"

template< typename T1_, typename T2_ >
inline typename element_type_of< T1_ >::type distance(
        T1_ const& p1, T2_ const p2,
        typename boost::enable_if<
            typename boost::mpl::and_<
                is_vector3<T1_>,
                is_vector3<T2_> >::type>::type* = 0)
{
    return std::sqrt(
        gsl_pow_2( p1[0] - p2[0] )
        + gsl_pow_2( p1[1] - p2[1] ) 
        + gsl_pow_2( p1[2] - p2[2] ) );
}

template<typename T_>
inline typename element_type_of<T_>::type distance(T_ const& p1, T_ const& p2)
{
    return distance(p1, p2, (void*)0);
}

template<typename T_>
inline T_ normalize(T_ const& p, 
                     typename element_type_of< T_ >::type const& r)
{
    typedef typename element_type_of< T_ >::type element_type;

    return multiply(p, r / length(p));
}


/**
 * Transpose the position pos1 so that it can be used with another 
 * position pos2.
 *
 * pos1 is transposed into one of mirror images of the cyclic boundary
 * condition so that the distance between pos1 and pos2 is smallest.
 *
 * Both of given pos1 and pos2 must be within the cyclic boundary.  However,
 * note that the returned transposed pos1 may not be within the cyclic boundary.
 */
template<typename T_>
inline T_ cyclic_transpose(T_ const& p0, T_ const& p1, T_ const& world_size, typename boost::enable_if<is_scalar<T_> >::type*)
{
    const T_ diff(p1 - p0), half(world_size / 2);
    if (diff > half)
    {
        return p0 + world_size;
    }
    else if (diff < -half)
    {
        return p0 - world_size;
    }
    else
    {
        return p0;
    }
}

template<typename T_>
inline T_ cyclic_transpose(T_ const& p0, T_ const& p1, typename element_type_of<T_>::type const& world_size, typename boost::enable_if<is_vector3<T_> >::type*)
{
    T_ retval;
    retval[0] = cyclic_transpose(p0[0], p1[0], world_size, (void*)0);
    retval[1] = cyclic_transpose(p0[1], p1[1], world_size, (void*)0);
    retval[2] = cyclic_transpose(p0[2], p1[2], world_size, (void*)0);
    return retval;
}

template<typename T1_, typename T2_>
inline T1_ cyclic_transpose(T1_ const& p0, T1_ const& p1, T2_ const& world_size)
{
    return cyclic_transpose(p0, p1, world_size, (void*)0);
}

template<typename T_>
inline T_ apply_boundary(T_ const& p1, T_ const& world_size, typename boost::enable_if<is_scalar<T_> >::type*)
{
    return modulo(p1, world_size);
}

template<typename T_>
inline T_ apply_boundary(T_ const& p1, 
                         typename element_type_of<T_>::type const& world_size,
                         typename boost::enable_if<is_vector3<T_> >::type*)
{
    return modulo(p1, world_size);
}

template<typename T1_, typename T2_>
inline T1_ apply_boundary(T1_ const& p1, T2_ const& world_size)
{
    return apply_boundary(p1, world_size, (void*)0);
}

#endif /* GEOMETRY_HPP */

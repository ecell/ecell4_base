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
                is_vector3<T2_> > >::type* = 0)
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
inline T_ normalize(T_ const& p)
{
    return divide(p, length(p));
}

template<typename T_>
inline T_ normalize(T_ const& p, 
                     typename element_type_of< T_ >::type const& r)
{
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

template<typename T1_, typename T2_>
inline typename element_type_of<T1_>::type distance_cyclic(
        T1_ const& p1, T2_ const& p2,
        typename element_type_of<T1_>::type const& world_size,
        typename boost::enable_if<
            typename boost::mpl::and_<
                is_vector3<T1_>,
                is_vector3<T2_> > >::type* = 0)
{
    return distance(p1, cyclic_transpose(p2, p1, world_size));
}

template<typename T_>
inline typename element_type_of<T_>::type
distance_cyclic(T_ const& p1, T_ const& p2,
                typename element_type_of<T_>::type const& world_size)
{
    return distance_cyclic(p1, p2, world_size, (void*)0);
}

template<typename T>
inline T spherical_to_cartesian(T const& s)
{
    typename element_type_of<T>::type const sintheta(std::sin(s[1]));
    T retval;
    retval[0] = s[0] * std::cos(s[2]) * sintheta;
    retval[1] = s[0] * std::sin(s[2]) * sintheta;
    retval[2] = s[0] * std::cos(s[1]);
    return retval;
}

template<typename T1, typename T2>
inline T1 rotate_vector(T1 const& v, T2 const& axis, double angle)
{
    double const c(std::cos(angle)), s(std::sin(angle)), cc(1. - c);
    double const mat[3][3] = {
        {
            c + cc * axis[0] * axis[0],
            cc * axis[0] * axis[1] - axis[2] * s,
            cc * axis[0] * axis[2] + axis[1] * s
        },
        {
            cc * axis[0] * axis[1] + axis[2] * s,
            c + cc * axis[1] * axis[1],
            cc * axis[1] * axis[2] - axis[0] * s
        },
        {
            cc * axis[0] * axis[2] - axis[1] * s,
            cc * axis[1] * axis[2] + axis[0] * s,
            c + cc * axis[2] * axis[2]
        }
    };

    return multiply(mat, v);
}

#endif /* GEOMETRY_HPP */

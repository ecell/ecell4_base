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
        pow_2( p1[0] - p2[0] )
        + pow_2( p1[1] - p2[1] ) 
        + pow_2( p1[2] - p2[2] ) );
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
inline T_ periodic_transpose(T_ const& p0, T_ const& p1, T_ const& world_size, typename boost::enable_if<is_scalar<T_> >::type*)
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
inline T_ periodic_transpose(T_ const& p0, T_ const& p1, typename element_type_of<T_>::type const& world_size, typename boost::enable_if<is_vector3<T_> >::type*)
{
    T_ retval;
    retval[0] = periodic_transpose(p0[0], p1[0], world_size, (void*)0);
    retval[1] = periodic_transpose(p0[1], p1[1], world_size, (void*)0);
    retval[2] = periodic_transpose(p0[2], p1[2], world_size, (void*)0);
    return retval;
}

template<typename T_>
inline T_ periodic_transpose(T_ const& p0, T_ const& p1, T_ const& edge_lengths, typename boost::enable_if<is_vector3<T_> >::type*)
{
    T_ retval;
    retval[0] = periodic_transpose(p0[0], p1[0], edge_lengths[0], (void*)0);
    retval[1] = periodic_transpose(p0[1], p1[1], edge_lengths[1], (void*)0);
    retval[2] = periodic_transpose(p0[2], p1[2], edge_lengths[2], (void*)0);
    return retval;
}

template<typename T1_, typename T2_>
inline T1_ periodic_transpose(T1_ const& p0, T1_ const& p1, T2_ const& world_size)
{
    return periodic_transpose(p0, p1, world_size, (void*)0);
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
inline T1_ apply_boundary(T1_ const& p1, T2_ const& edge_lengths, typename boost::enable_if<typename boost::mpl::and_<is_vector3<T1_>, is_vector3<T2_> > >::type*)
{
    return modulo(p1, edge_lengths);
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
    return distance(p1, periodic_transpose(p2, p1, world_size));
}

template<typename T1_, typename T2_, typename T3_>
inline typename element_type_of<T1_>::type distance_cyclic(
        T1_ const& p1, T2_ const& p2, T3_ const& edge_lengths,
        typename boost::enable_if<
            typename boost::mpl::and_<
                is_vector3<T1_>,
                is_vector3<T2_>,
                is_vector3<T3_> > >::type* = 0)
{
    return distance(p1, periodic_transpose(p2, p1, edge_lengths));
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

// reflect 
template<typename coordT>
coordT reflect_plane(const coordT& begin, const coordT& end,
                     const coordT& normal, const coordT& plane)
{
    typedef typename element_type_of<coordT>::type valueT;
//     assert(std::abs(length(normal) - 1.0) < 1e-12);
    const valueT norm_b = dot_product((begin - plane), normal);
    const valueT norm_e = dot_product((end - plane), normal);
    if(norm_b == 0.0)
    {
        throw std::invalid_argument("reflection: begin is on the plane");
    }
    else if(norm_b * norm_e > 0.0 && std::abs(norm_e) < 1e-10)
    {
        return (begin * 1e-10) + (end * (1.0 - 1e-10));
    }
    else if(norm_b * norm_e < 0.0 && std::abs(norm_e) < 1e-10)
    {
        return begin * 1e-10 + (end - (normal * (norm_e * 2.0))) * (1. - 1e-10);
    }
    else if(norm_b * norm_e > 0.0)
    {
        return end;
    }
    else
    {
        return end - (normal * (norm_e * 2.0));
    }
}

template<typename coordT>
inline typename element_type_of<coordT>::type
angle(const coordT& lhs, const coordT& rhs)
{
    typedef typename element_type_of<coordT>::type valueT;
    const valueT lensq_l = length_sq(lhs);
    const valueT lensq_r = length_sq(rhs);
    const valueT inner = dot_product(lhs, rhs);
    return acos(inner / std::sqrt(lensq_l * lensq_r));
}


#endif /* GEOMETRY_HPP */

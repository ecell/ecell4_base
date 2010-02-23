#ifndef LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_HPP

#include <algorithm>
#include <cmath>
#include <gsl/gsl_pow_int.h>
#include <boost/mpl/bool.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/array.hpp>

#include "utils/array_traits.hpp"

#define CREATE_VECTOR_LIMIT_REPEAT 16

template<typename T_, std::size_t N_>
struct is_vector: public boost::mpl::false_ {};

template<typename T_, std::size_t N_>
struct is_vector<boost::array<T_, N_>, N_>: public boost::mpl::true_ {};

template<typename T_>
struct is_scalar: public boost::is_arithmetic<T_> {};

template<typename T_>
struct is_vector2: public is_vector<T_, 2> {};

template<typename T_>
struct is_vector3: public is_vector<T_, 3> {};


template<typename T_>
inline T_ add( T_ const& p1, T_ const& p2, typename boost::enable_if<is_scalar<T_> >::type* = 0)
{
    return p1 + p2;
}

template<typename T_>
inline T_ subtract( T_ const& p1, T_ const& p2, typename boost::enable_if<is_scalar<T_> >::type* = 0)
{
    return p1 - p2;
}

template<typename T_>
inline T_ multiply( T_ const& p1, T_ const& p2, typename boost::enable_if<is_scalar<T_> >::type* = 0)
{
    return p1 * p2;
}

template<typename T_>
inline T_ divide( T_ const& p1, T_ const& p2, typename boost::enable_if<is_scalar<T_> >::type* = 0)
{
    return p1 / p2;
}

template<typename T_>
inline T_ modulo( T_ const& p1, T_ const& p2 )
{
    T_ r = p1 % p2;
    if (r != 0 && (r > 0) == (p2 < 0))
        r += p2;
    return r;
}

template<>
inline float modulo( float const& p1, float const& p2 )
{
    float r = std::fmod(p1, p2);
    if (r != 0 && (r > 0) == (p2 < 0))
        r += p2;
    return r;
}

template<>
inline double modulo( double const& p1, double const& p2 )
{
    double r = std::fmod(p1, p2);
    if (r != 0 && (r > 0) == (p2 < 0))
        r += p2;
    return r;
}

template<typename T_>
inline T_ negate(T_ const& v, typename boost::enable_if<is_scalar<T_> >::type* = 0)
{
    return -v;
}

template<typename T_>
inline T_ abs(T_ const& v, typename boost::enable_if<is_scalar<T_> >::type* = 0)
{
    return std::fabs(v);
}

template< typename T_ >
inline T_ add( T_ const& p1, T_ const& p2, typename boost::enable_if<is_vector3<T_> >::type* = 0 )
{
    T_ retval;
    retval[0] = add(p1[0], p2[0]);
    retval[1] = add(p1[1], p2[1]);
    retval[2] = add(p1[2], p2[2]);
    return retval;
}

template< typename T_ >
inline T_ subtract( T_ const& p1, T_ const& p2, typename boost::enable_if<is_vector3<T_> >::type* = 0 )
{
    T_ retval;
    retval[0] = subtract(p1[0], p2[0]);
    retval[1] = subtract(p1[1], p2[1]);
    retval[2] = subtract(p1[2], p2[2]);
    return retval;
}

template<typename T_>
inline T_ divide( T_ const& p1, typename T_::value_type p2, typename boost::enable_if<is_vector3<T_> >::type* = 0 )
{
    T_ retval;
    retval[0] = divide(p1[0], p2);
    retval[1] = divide(p1[1], p2);
    retval[2] = divide(p1[2], p2);
    return retval;
}

template<typename T_>
inline T_ multiply( T_ const& p1, typename T_::value_type p2, typename boost::enable_if<is_vector3<T_> >::type* = 0 )
{
    T_ retval;
    retval[0] = multiply(p1[0], p2);
    retval[1] = multiply(p1[1], p2);
    retval[2] = multiply(p1[2], p2);
    return retval;
}

template<typename T_>
inline T_ modulo( T_ const& p1, typename T_::value_type p2, typename boost::enable_if<is_vector3<T_> >::type* = 0 )
{
    T_ retval;
    retval[0] = modulo(p1[0], p2);
    retval[1] = modulo(p1[1], p2);
    retval[2] = modulo(p1[2], p2);
    return retval;
}

template<typename T_>
inline T_ negate(T_ const& v, typename boost::enable_if<is_vector3<T_> >::type* = 0)
{
    T_ retval;
    retval[0] = -v[0];
    retval[1] = -v[1];
    retval[2] = -v[2];
    return retval;
}

template<typename T_>
inline T_ abs(T_ const& v, typename boost::enable_if<is_vector3<T_> >::type* = 0)
{
    T_ retval;
    retval[0] = abs(v[0]);
    retval[1] = abs(v[1]);
    retval[2] = abs(v[2]);
    return retval;
}

template<typename T_>
inline typename element_type_of<T_>::type dot_product(T_ const& p1, T_ const& p2, typename boost::enable_if<is_vector3<T_> >::type* = 0)
{
    return multiply(p1[0], p2[0])
           + multiply(p1[1], p2[1])
           + multiply(p1[2], p2[2]);
}

template<typename T_>
inline typename element_type_of<T_>::type length_sq(T_ const& r, typename boost::enable_if<is_vector2<T_> >::type* = 0)
{
    return gsl_pow_2(r[0]) + gsl_pow_2(r[1]);
}

template< typename T_ >
inline typename element_type_of< T_ >::type length_sq(T_ const& r, typename boost::enable_if<is_vector3<T_> >::type* = 0)
{
    return gsl_pow_2(r[0]) + gsl_pow_2(r[1]) + gsl_pow_2(r[2]);
}

template< typename T_ >
inline typename element_type_of< T_ >::type length(T_ const& r)
{
    return std::sqrt(length_sq(r));
}

#define CREATE_VECTOR_INNER_TPL(__z__, __n__, __d__) \
    __d__[__n__] = BOOST_PP_CAT(p, __n__);

#define CREATE_VECTOR_TPL(__z__, __n__, __d__) \
template<typename T_> \
inline T_ create_vector(\
        BOOST_PP_ENUM_PARAMS(__n__, typename element_type_of<T_>::type const& p), \
        typename boost::enable_if<is_vector<T_, __n__> >::type* = 0) \
{ \
    T_ retval; \
    BOOST_PP_REPEAT_ ## __z__(__n__, CREATE_VECTOR_INNER_TPL, retval) \
    return retval; \
}

BOOST_PP_REPEAT_FROM_TO(2, CREATE_VECTOR_LIMIT_REPEAT, CREATE_VECTOR_TPL, )

#undef CREATE_VECTOR_TPL
#undef CREATE_VECTOR_INNER_TPL

#endif /* LINEAR_ALGEBRA_HPP */
